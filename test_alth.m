             
%% step 2: design branch-and-bound method.
    tau_init=[0,0,1,0,0,0];
    level=3;
    gap=5;  
    max_rotation=args.branch_max_rotation;  
    max_skew=args.branch_max_skew;
    candidate_matrix=cell(3, 2*args.branch_accuracy+1);    
    BLACK_MATRIX=zeros(focus_size(1)*level+gap*(level-1), focus_size(2)*(2*args.branch_accuracy+1)+gap*2*args.branch_accuracy);
    %% step 3: begin branch-and-bound
    normal_outer_max_iter=args.outer_max_iter;
    normal_display_inter=args.display_result;
    args.outer_max_iter=1; % for debug, set it to 1;
    args.display_result=0;  
  
    for i=1:level
        for j=1:2*args.branch_accuracy+1
            tau=tau_init;
            candidate_matrix{1, j}=eye(3);
            theta=-max_rotation+(i-1)*max_rotation/args.branch_accuracy;
            tau(i+3)=theta;
            candidate_matrix{i,j}=tau;
           args.figure_no=(i-1)*level+j;
           args.save_path=[];
%            [Dotau, A, E, f, tfm_matrix, error_sign, UData, VData, XData, YData, A_scale]=tilt_kernel(input_image, args.mode, center, focus_size, tfm_matrix, args);
           %% update 1: assume that there's no great error, we directly
           %% compute nuclear norm of Dotau as the selection criterion;
           %% record result
            image_size=size(input_image);
            image_center=floor(center);
            focus_center=zeros(2, 1);
            focus_center(1)=floor((1+focus_size(2))/2);
            focus_center(2)=floor((1+focus_size(1))/2);
            UData=[1-image_center(1) image_size(2)-image_center(1)];
            VData=[1-image_center(2) image_size(1)-image_center(2)];
            XData=[1-focus_center(1) focus_size(2)-focus_center(1)];
            YData=[1-focus_center(2) focus_size(1)-focus_center(2)];
            Dotau= generate_Branch_and_Bound(input_image,tau,XData,YData,image_center);
%             tfm=fliptform(maketform('projective', tfm_matrix'));
%             Dotau=imtransform(input_image, tfm, 'bilinear', 'XData', XData, 'YData', YData, 'UData', UData, 'VData', VData, 'Size', focus_size);
            Dotau=Dotau/norm(Dotau, 'fro');
            [U S V]=svd(Dotau);
            f=sum(sum(S));
           disp(['branching: level=', num2str(i), ', idx=', num2str(j)]);
           start=[(focus_size(1)+gap)*(i-1)+1, (focus_size(2)+gap)*(j-1)+1];
           BLACK_MATRIX(start(1):(start(1)+focus_size(1)-1), start(2):(start(2)+focus_size(2)-1))=Dotau;
           f_branch(i, j)=f;
%            A_branch{i, j}=A;
%            E_branch{i, j}=E;
           Dotau_branch{i, j}=Dotau;
           result_tau{i, j}=tau;
        end
        [value index]=min(f_branch(i, :));
        tau_init=result_tau{i, index};
    end

  
    %% step 4: show inter result if necessary.
    if normal_display_inter==1
        showimage(BLACK_MATRIX, 98);
    end
    args.outer_max_iter=normal_outer_max_iter;
    args.display_result=normal_display_inter;