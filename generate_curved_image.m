function curved_image = generate_curved_image( img, tau )

Img0 = im2double(img);

% Parameters
R = tau(1);
theta = tau(2);
h = tau(3);


[ H, W ] = size( Img0 );

% Refine Paramters

h0 = ( H - 1 )/2;

%

[ X0, Y0 ] = meshgrid( 1 : W, -h0 : h0 );

   

XI = ( X0 + (Y0-h)*tan(theta) ) - ( R - ( R^2 - ((Y0-h)/cos(theta)).^2 ).^0.5 );
YI = R*asin((Y0-h)/(R*cos(theta)));

Img_out2 = interp2( X0, Y0, Img0, XI, YI, 'cubic' );
Img_out2( isnan( Img_out2 ) ) = 1;

curved_image = Img_out2;