function [Dotau, du, dv] = generate_Dotau_curve2(input_image,tau,input_du,input_dv,UData, VData, XData,YData)




% Parameters
R = tau(1);
theta = tau(2);
h = tau(3);% ºŸ…Ë√ª”–∆´“∆

[ X0, Y0 ] = meshgrid( UData(1) : UData(2) , VData(1) : VData(2) );


%% Method II

XII = X0 + R*( 1 - cos( (Y0) / R ) ) - R*sin( (Y0)/R )*sin(theta);
YII = R*sin( (Y0)/R )*cos(theta)+h;

Img_out = interp2( X0, Y0, input_image, XII, YII, 'cubic' );
Img_out( isnan( Img_out ) ) = 1;

du0 = interp2( X0, Y0, input_du, XII, YII, 'cubic' );
du0( isnan( du0 ) ) = 1;

dv0 = interp2( X0, Y0, input_dv, XII, YII, 'cubic' );
dv0( isnan( dv0 ) ) = 1;

Dotau = Img_out( YData(1) : YData(2), XData(1) : XData(2) );
du = du0( YData(1) : YData(2), XData(1) : XData(2) );
dv = dv0( YData(1) : YData(2), XData(1) : XData(2) );