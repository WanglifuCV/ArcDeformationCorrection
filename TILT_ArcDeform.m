function [Dotau, A, E, f,tau, focus_size, error_sign,A_scale]=TILT_ArcDeform(varargin)

% TILT will align the drawn-part of the image to its frontal.
% TILT is built on the kernel_component tilt_kernel.
% TILT(input_image, mode, based_points, focus_size) is the simplest form of
% the parameter, and these four parameters must be specified. There are
% also many other optional parameters.
%
% ----------------------necessary input------------------------------------
% input_image:  height-by-width real matrix or height*width*3 real matrix but
%               we will only preserve the first channel for the second
%               case.
% mode:         one of 'euclidean', 'affine', 'homography'.
%
% [note]: One of the following sets of parameters must be specified.
% Set 1: initial_points:
%               2-by-1 real vector in x-y co-ordinates, top-left point of
%               the focus. In this situation our algorithm will
%               automatically decide the UData, VData, XData and YData as
%               the output.
% Set 2: UData, VData, XData and YData.
%               All 1-by-2 real vector. These specifying the coordinate of
%               the input image and the transformed image, and we'll work
%               in this coordinate system by default.
%
% ----------------------optional parameters--------------------------------
% INITIAL_MATRIX:   3-by-3 initialization matrix for tfm_matrix. If this is
%                   not set, tfm_matrix will be initialized to identity.
% BRANCH:           0 or 1, 1 for turnning on branch-and-bound, but this is
%                   only allowed in the AFFINE case.
% BLUR:             0 or 1, 1 for turnning on BLUR.
% PYRAMID:          0 or 1, 1 for turnning on pyramid.
% NO_TRANSLATION:   0 or 1, 1 for no translation.
%
% INNER_TOL:        positive real value, inner_loop threshold.
% INNER_C:          positive real value, inner_loop lambda=C/sqrt(m).
% INNER_MU:         positive real value, inner_loop for ALM mu.
% INNER_MAXITER:    positive integer, maximum iteration for inner_loop.
% INNER_DISPLAY:    positive integer, whether we display the results.
% OUTER_TOL:        positive real value, outer_loop threshold.
% OUTER_MAXITER:    positive integer, maximum iteration for ouer_loop.
% OUTER_DISPLAY:    positive integer, whether we display the results for
%                   the outer loop.
% FOCUS_THRE:       positive integer, smallest edge length threshold in
%                   pyramid
% OUTER_TOL_STEP:   positive real value, We relax threshold for
%                   outer-loop each time we move downstairs in pyramid
%                   by tol=tol*OUTER_TOL_STEP
% BLUR_SIGMA:       positive real value, standard derivation of the blur
%                   kernel.
% BLUR_NEIGHBOR:    positive integer, size of effective blur neighbourhood.
% BRANCH_MAXITER:   positive integer, we need have extremely high accuracy
%                   in branch-and=bound. So we separately set
%                   branch_maxiter for branch_and_bound.
% BRANCH_ACCURACY:  positive integer, we split the whole parameter region
%                   to search into 2*R_accuracy+1 sub-region.
% BRANCH_MAX_ROTATION:
%                   positive real, by default pi/8, specifying how large to
%                   do branch-and-bound for rotation.
% BRANCH_MAX_SKEW:
%                   positive, real, by default
% PYRAMID_MAXLEVEL: positive integer, we only run TILT on the highest
%                   PYRAMID_MAXLEVEL levels in the pyramid.
% DISPLAY_INTER:    0 or 1, whether we display results for every call to
%                   tilt_kernel.m
% FOCUS_SIZE:       1*2 row vector, forcing focus_size to be something.
% SAVE_PATH:        string, full path of where to save the results.
%
% -------------------------output------------------------------------------
% Dotau:        real matrix with same size as focus_size, aligned images.
% A:            low-rank part of Dotau;
% E:            sparse-error part of Dotau;
% f:            value of objective-function;
% tfm_matrix:   resulted transform matrix.
% focus_size:   1*2 positive integer vector, size of the focus window in
%               r-c coordinate.
% error_sign:   0 or 1, 1 for trival solutions.
% UData, VData: 1*2 real vector, position of the input image.
% XData, YData: 1*2 real vector, position of the transformed image.


tic;
args=parse_inputs(varargin{:});
%% Automatically resize and cut the original image to avoid unnecessary
%% computation that happens far away from the input image.
%%% back up
original_args=args;


%%% cut the image
expand_rate=0.8;
initial_points=args.initial_points;   %左上、右下
left_bound=ceil(max(initial_points(1, 1)-expand_rate*(initial_points(1, 2)-initial_points(1, 1)), 1));  %1
right_bound=floor(min(initial_points(1, 2)+expand_rate*(initial_points(1, 2)-initial_points(1, 1)), size(args.input_image, 2))) ; %1280
top_bound=ceil(max(initial_points(2, 1)-expand_rate*(initial_points(2, 2)-initial_points(2, 1)), 1)) ;                          %1
bottom_bound=floor(min(initial_points(2, 2)+expand_rate*(initial_points(2, 2)-initial_points(2, 1)), size(args.input_image, 1))) ; %1024
new_image=zeros(bottom_bound-top_bound+1, right_bound-left_bound+1,  size(args.input_image, 3));

for c=1:size(args.input_image, 3)
    new_image(:, :, c)=args.input_image(top_bound:bottom_bound, left_bound:right_bound, c);
end
args.input_image=uint8(new_image);
args.center=args.center+[1-left_bound; 1-top_bound];  %center of the cut image




%%% down-sample the image if the focus is too large.
pre_scale_matrix=eye(3);
focus_threshold=800; % The max length of the focus is 100. if too large ,downsample to matrix of this size.

min_length=min(original_args.focus_size);
if min_length>focus_threshold
    s=min_length/focus_threshold;
    dst_pt_size=round([bottom_bound-top_bound+1, right_bound-left_bound+1]/s);
    % %     args.input_image=imfilter(args.input_image, fspecial('gaussian', ceil(2^(s-1)), ceil(2^(s-1))));
    args.input_image=imresize(args.input_image, dst_pt_size);
    
    pre_scale_matrix=diag([s s 1]);
end


pre_scale=pre_scale_matrix(1,1);
%%% adjust the parsed parameters.
initial_tfm_matrix=args.initial_tfm_matrix;
initial_tfm_matrix=inv(pre_scale_matrix)*initial_tfm_matrix*pre_scale_matrix;
args.initial_tfm_matrix=initial_tfm_matrix;

args.focus_size=round(args.focus_size/pre_scale_matrix(1, 1));
args.center=args.center/pre_scale_matrix(1, 1);
parent_path=args.save_path;

if ~exist(args.save_path) || ~isdir(args.save_path)
    mkdir(args.save_path);
end

img=args.input_image;
img = double(adapthisteq(uint8(img)));
if args.branch==1
    %% Branch-and-bound if set
    %% step 1: prepare data for the lowest resolution.
    downsample_rate=0.5;
    total_scale=floor(log2(min(args.focus_size)/args.focus_threshold));
    downsample_matrix=[0.5 0 0; 0 0.5 0; 0 0 1];
    scale_matrix=downsample_matrix^total_scale;
    branch_scale=0.5^total_scale*pre_scale;
    tfm=maketform('projective', scale_matrix');
    if args.blur
        %% blur if set
        input_image=imfilter(args.input_image, fspecial('gaussian', ceil(args.blur_kernel_size_k*2^total_scale), ceil(args.blur_kernel_sigma_k*2^total_scale)));
    else
        input_image=args.input_image;
    end
    input_image=imtransform(input_image, tfm, 'bicubic');
    if size(input_image, 3)>1
        input_image=rgb2gray(input_image);
    end
    input_image=double(input_image);
    % initial_tfm_matrix=scale_matrix*args.initial_tfm_matrix*inv(scale_matrix);
    center=floor(transform_point(args.center, scale_matrix));
    
    focus_size=floor(args.focus_size/2^total_scale);
    
    input_image_threshold = im2bw( input_image );
    
    [ y_index, ~ ] = find(input_image_threshold == 0 );
    
    y_center = mean( y_index );
    
%     figure(50);imshow( input_image, [] );
%     hold on;plot( [ 1  size( input_image, 2 ) ], [ y_center y_center ], 'r' );
    
    %% step 2: design branch-and-bound method.
    
    level=(args.theta_max/args.theta_accuracy);
    gap=5;
    
    f_branch=zeros(2*level+1, args.R_accuracy);
    Dotau_branch=cell(2*level+1, args.R_accuracy);
    result_tau=cell(2*level+1, args.R_accuracy);
    
    R_min = 1.1*size(input_image,1)/2;
    R_max = 10*size(input_image,1)/2;
    R_gap = (R_max/R_min).^(1/args.R_accuracy);
    
    tau_init=[0 0 0];
    
%     BLACK_MATRIX=zeros(focus_size(1)*level+gap*(level-1), focus_size(2)*(args.R_accuracy)+gap*2*args.R_accuracy);
    %% step 3: begin branch-and-bound
    normal_outer_max_iter=args.outer_max_iter;
    normal_display_inter=args.display_result;
    args.outer_max_iter=1; % for debug, set it to 1;
    args.display_result=0;
    
    tau=tau_init;
    
    for i= -level : level
        tau(2) = args.theta_accuracy * i;
        for j=1:args.R_accuracy
            
            R = R_min*R_gap^(j-1);
            tau(1) = R;
            
%             args.figure_no=(i-1)*level+j;
            args.save_path=[];
            %            [Dotau, A, E, f, tfm_matrix, error_sign, UData, VData, XData, YData, A_scale]=tilt_kernel(input_image, args.mode, center, focus_size, tfm_matrix, args);
            %% update 1: assume that there's no great error, we directly
            %% compute nuclear norm of Dotau as the selection criterion;
            %% record result
            image_size=size(input_image);
            image_center=floor(center);
            focus_center=zeros(2, 1);
            
            UData=[ 1 image_size(2) ];
            VData=[ -image_center(2) image_size(1)-image_center(2)- 1];
            XData=[ image_center(1) - floor(focus_size(2)/2)+1 image_center(1) + ( focus_size(2)-floor( focus_size(2)/2 ) ) ];
            YData=[ image_center(2) - floor(focus_size(1)/2)+1 image_center(2) + ( focus_size(1)-floor( focus_size(1)/2 ) ) ];
            
            Dotau= generate_Branch_and_Bound2(input_image,tau, UData, VData, XData, YData );
            
            
            figure(350);imshow(uint8(Dotau));title( num2str( tau(1) ) );
            Dotau=Dotau/norm(Dotau, 'fro');
            [~, S, ~]=svd(Dotau);
            f=sum(sum(S));
            disp(['branching: level=', num2str(i+level+1), ', idx=', num2str(j)]);
%             start=[(focus_size(1)+gap)*(i-1)+1, (focus_size(2)+gap)*(j-1)+1];
%             BLACK_MATRIX(start(1):(start(1)+focus_size(1)-1), start(2):(start(2)+focus_size(2)-1))=Dotau;
            f_branch(i+level+1, j)=f;
            Dotau_branch{i+level+1, j}=Dotau;
            result_tau{i+level+1, j}=tau;
        end
        %         [~, index]=min(f_branch(i, :));
        %         tau_init=result_tau{i, index};
        
    end
    [i_index, j_index] = find( f_branch == min( f_branch(:) ) );
    tau_init = result_tau{ i_index, j_index };
    
    Dotau_init = generate_Branch_and_Bound2(input_image,tau_init, UData, VData, XData, YData );
    figure(110);imshow( Dotau_init,[] );title('branch-and-bound');
    
    %% step 4: show inter result if necessary.
    if normal_display_inter==1
        %    showimage(BLACK_MATRIX, 98);
    end
    args.outer_max_iter=normal_outer_max_iter;
    args.display_result=normal_display_inter;
end

tau_init(1) = tau_init(1)/(pre_scale*downsample_rate^total_scale);

%% Do pyramid if necessary
if args.pyramid==1
    %% define parameters
    downsample_rate=0.5;
    % upsample_rate=2;
    downsample_matrix=[0.5 0 0; 0 0.5 0; 0 0 1];
    %upsample_matrix=inv(downsample_matrix);
    %     total_scale=ceil(max(log2(min(args.focus_size)/args.focus_threshold), 0))
    for scale=total_scale:-1:0
        downsample_scale = downsample_rate^scale;
        %upsample_scale=upsample_rate^scale;
        %% begin each level of the pyramid
        if total_scale-scale>=args.pyramid_max_level
            break;
        end
        %% Blur if required
        if args.blur==1 && scale~=0
            input_image=imfilter(img, fspecial('gaussian', ceil(args.blur_kernel_size_k*2^scale), ceil(args.blur_kernel_sigma_k*2^scale)));
        else
            input_image=args.input_image;
        end
        
        %% prepare image and initial tfm_matrix
        img_scale = pre_scale* downsample_scale;
        scale_matrix=downsample_matrix^scale;
        tfm=maketform('projective', scale_matrix');
        input_image=imtransform(input_image, tfm, 'bicubic');
        %size(input_image)
        %figure(11);imshow(uint8(input_image));
        %        tfm_matrix=scale_matrix*initial_tfm_matrix*inv(scale_matrix);
        center=floor(transform_point(args.center, scale_matrix));
        focus_size=floor(args.focus_size/2^scale);
        args.save_path=fullfile(parent_path, ['pyramid', num2str(scale)]);
        args.figure_no=100+total_scale-scale+1;
        tau(1)=tau_init(1)*(pre_scale*downsample_rate^scale);
        tau(2)=tau_init(2);
        tau(3)=tau_init(3)*(pre_scale*downsample_rate^scale);
        [Dotau, A, E, f,tau, error_sign]=tilt_kernel(input_image, center, focus_size,tau,img_scale,args);
        %% update tfm_matrix of the highest-resolution level.
        % initial_tfm_matrix=inv(scale_matrix)*tfm_matrix*scale_matrix;
        
        tau_init(1)=tau(1)/(pre_scale*downsample_rate^scale);
        tau_init(2)=tau(2);
        tau_init(3)=tau(3)/(pre_scale*downsample_rate^scale);
        args.outer_tol=args.outer_tol*args.outer_tol_step;
    end
    % tfm_matrix=initial_tfm_matrix;
    %     Rectify_image=generate_new_image_curve2(input_image,tau,floor(center),scale);

    tau=tau_init;
else
    %% No Pyramid
    %      tau_init = [0; 0; 1; 0; 0; 0];
    tau=tau_init;
    %% Blur if required
    if args.blur==1
        img_size=size(args.input_image);
        img_size=img_size(1:2);
        
        input_image=imfilter(args.input_image, fspecial('gaussian', ceil(args.blur_kernel_size_k*max(img_size)/50), ceil(args.blur_kernel_sigma_k*max(img_size)/50)));
        %input_image=imfilter(img, fspecial('gaussian', ceil(args.blur_kernel_size_k*max(img_size)/50), ceil(args.blur_kernel_sigma_k*max(img_size)/50)));
    else
        input_image=args.input_image;
    end
    args.figure_no=101;
    args.save_path=fullfile(parent_path, 'pyramid0');
    [Dotau, A, E, f,tau,error_sign]=tilt_kernel(input_image, args.center, args.focus_size,tau, 1,args);
end
toc

%%
%   figure; imshow(A, [], 'DisplayRange', [0 max(max(max(A)))]); title('Low^R_ank');


A_scale=pre_scale;
% scale_matrix=[1 0 0 ;0 1 0;0 0 1];
% image_center=floor(transform_point(args.center, scale_matrix));



function args=parse_inputs(varargin)
iptchecknargin(3,Inf,nargin,mfilename);
%% default value
args.input_image=varargin{1};
args.mode=varargin{2};
args.initial_points=varargin{3};
args.initial_tfm_matrix=eye(3);
args.outer_tol=1e-4;
args.outer_max_iter=50;
args.outer_display_period=1;
args.inner_tol=1e-4;
args.inner_c=1;
args.inner_mu=[];
args.inner_display_period=100;
args.inner_max_iter=inf;
args.blur=0;

args.pyramid=1;
args.branch=1;

args.focus_threshold=50; %原来是50 % when doing pyramid the smallest focus_edge we can tolerate.
args.outer_tol_step=10; % as resolution goes high how relaxed should the outer_tol be.
args.blur_kernel_size_k=2;% neighbourhood scalar for the size of the blur kernel.
args.blur_kernel_sigma_k=2;% standard derivation scalar for blur kernel.
args.pyramid_max_level=4; % number of pyramid levels we want to act on.
args.branch_max_iter=10; % in each branch, how much iteration we take.
args.R_accuracy=5; % higher means smaller step-width.
args.theta_accuracy = 2*pi/180;% 以2°为间隔
args.theta_max = 20*pi/180;%最大20°
args.display_result=0;
args.focus_size=[];
args.save_path=[];

args.R_accuracy = 21;

if strcmpi(args.mode, 'affine')
    args.branch_max_rotation=pi/9;
    args.branch_max_skew=0.7; % purely empirical, just to emphasize that in affine mode, the skew should be smaller...
else
    args.branch_max_rotation=pi/6;
    args.branch_max_skew=1;
end
for i=4:2:nargin
    switch(lower(varargin{i}))
        case 'outer_tol'
            args.outer_tol=varargin{i+1};
        case 'outer_maxiter'
            args.outer_max_iter=varargin{i+1};
        case 'outer_display'
            args.outer_display_period=varargin{i+1};
        case 'inner_tol'
            args.inner_tol=varargin{i+1};
        case 'inner_mu'
            args.inner_mu=varargin{i+1};
        case 'inner_display'
            args.inner_display_period=varargin{i+1};
        case 'inner_maxiter'
            args.inner_max_iter=varargin{i+1};
        case 'inner_c'
            args.inner_c=varargin{i+1};
        case 'blur'
            args.blur=varargin{i+1};
        case 'pyramid'
            args.pyramid=varargin{i+1};
        case 'branch'
            args.branch=varargin{i+1};
        case 'focus_thre'
            args.focus_threshold=varargin{i+1};
        case 'outer_tol_step'
            args.outer_tol_step=varargin{i+1};
        case 'blur_sigma'
            args.blur_kernel_sigma_k=varargin{i+1};
        case 'blur_neighbor'
            args.blur_kernel_size_k=varargin{i+1};
        case 'pyramid_maxlevel'
            args.pyramid_max_level=varargin{i+1};
        case 'branch_maxiter'
            args.branch_max_iter=varargin{i+1};
        case 'R_accuracy'
            args.R_accuracy=varargin{i+1};
        case 'no_translation'
            args.no_translation=varargin{i+1};
        case 'initial_matrix'
            args.initial_tfm_matrix=varargin{i+1};
        case 'display_inter'
            args.display_result=varargin{i+1};
        case 'initial_points'
            args.initial_points=varargin{i+1};
        case 'focus_size'
            args.focus_size=floor(varargin{i+1});
        case 'save_path'
            args.save_path=varargin{i+1};
        case 'branch_max_rotation'
            args.branch_max_rotation=varargin{i+1};
        case 'branch_max_skew'
            args.branch_max_skew=varargin{i+1};
    end
end

%% do some initialization
if size(args.initial_points, 2)==2
    args.initial_points=floor(args.initial_points);
    args.focus_size=[args.initial_points(2, 2)-args.initial_points(2, 1)+1 args.initial_points(1, 2)-args.initial_points(1, 1)+1];
    args.center=floor(mean(args.initial_points, 2));%对inital_points按行取平均，在inital_points中，左上角的点坐标存第一列，右下角点存第二列
    %截取的窗口的中心的坐标
elseif size(args.initial_points, 2)==4
    args.center=floor(mean(args.initial_points, 2));
    pt_mean=args.initial_points-args.center*ones(1, 4);
    if isempty(args.focus_size)
        args.focus_size=mean(abs(pt_mean), 2)*2;
        args.focus_size=floor([args.focus_size(2) args.focus_size(1)]);
    end
    focus_center=floor([1+args.focus_size(2) 1+args.focus_size(1)]/2);
    XData=[1-focus_center(1) args.focus_size(2)-focus_center(1)];
    YData=[1-focus_center(2) args.focus_size(1)-focus_center(2)];
    X=[XData(1) XData(2) XData(2) XData(1);...
        YData(1) YData(1) YData(2) YData(2)];
    args.initial_tfm_matrix=compute_homography(X, floor(pt_mean));
end