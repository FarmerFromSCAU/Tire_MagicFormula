%%
%Scaling Coefficients
%beta(1) LCX Scale factor of Fx shape factor
%beta(2) LMUX Scale factor of Fx peak friction coefficient
%beta(3) LEX Scale factor of Fx curvature factor
%beta(4) LKX Scale factor of slip stiffness
%beta(5) LHX Scale factor of Fxhorizontal shift
%beta(6) LVX Scale factor of Fx vertical shift
%beta(7) LXAL Scale factor of alpha influence on Fx
%%
%Longitudinal Coefficients
%beta(8) PCX1 Shape factor Cfx for longitudinal force
%beta(9) PDX1 Longitudinal friction Mux at Fznom
%beta(10) PDX2 Variation of friction Mux with load
%beta(11) PDX3 Variation of ftiction Mux with camber
%beta(12) PEX1 Longitudinal curvature Efx at Fznom
%beta(13) PEX2 Variation of curvature Efx with load
%beta(14) PEX3 Variation of curvature Efx with load squared
%beta(15) PEX4 Factor in curvature Efx while driving
%beta(16) PKX1 Longitudinal slip stiffness Kfx/Fz at Fznom
%beta(17) PKX2 Variation of slip stiffness Kfx/Fz with load
%beta(18) PKX3 Exponent in slip stiffness Kfx/Fz with load
%beta(19) PHX1 Horizontal shift Shx at Fznom
%beta(20) PHX2 Variation of shift Shx with load
%beta(21) PVX1 Vertical shift Svx/Fz at Fznom
%beta(22) PVX2 Variation of shift Svx/Fz with load
%beta(23) RBX1 Slope factor for combined slip Fx reduction
%beta(24) RBX2 Variation of slope Fx reduction with kappa
%beta(25) RBX3 Influence of camber on stiffness for Fx combined
%beta(26) RCX1 Shape factor for combined slip Fx reduction
%beta(27) REX1 Curvature factor of combined Fx
%beta(28) REX2 Curvature factor of combined Fx with load
%beta(29) RHX1 Shift factor for combined slip Fx reduction
%beta(30) PPX1 Linear pressure effect on slip stiffness
%beta(31) PPX2 Quadratic pressure effect on longitudinal friction
%beta(32) PPX3 Linear pressure effect on longitunal friction
%beta(33) PPX4 Quadratic pressure effect on longitudinal friction
%%
%TurnSlip Coefficients
%beta(34) PDXP1 Peak Fx reduction due to spin parameter
%beta(35) PDXP2 Peak Fx reduction due to spin wuth varying load parameter
%beta(36) PDXP3 Peak Fx reduction due to spin with kappa parameter
%%
%beta(37) LFZ0 Scale factor of nominal (rated) load
%beta(38) PECP1 Camber w.r.t spin reduction factor parameter in camber stiffness
%beta(39) PECP2 Camber w.r.t spin varying with load parameter in camber stiffness
%%
function Fx=LongitudinalForceFx(beta,x)
    Fz0=1112;
    Fz0a=-1*Fz0;
    dfz=(x(:,1)-Fz0a)/Fz0a;
    dpi=0;
    %Pure Slip
    Xi1=1;
    Shx=(beta(29)+beta(20)*dfz)*beta(5);
    KappaX=x(:,2)+Shx;
    Cx=beta(8)*beta(1);
    Mux=(beta(9)+beta(10)*dfz)*(1+beta(32)*dpi+beta(33)*dpi^2).*(1-beta(33).*x(:,3).^2)*beta(2);
    Ex=(beta(12)+beta(13)*dfz+beta(14)*dfz.^2).*(1-beta(15)*sign(KappaX))*beta(3);
    Kxk=x(:,1).*(beta(16)+beta(17)*dfz).*exp(beta(18)*dfz).*(1+beta(30)*dpi+beta(31)*dpi.^2)*beta(4);
    Dx=Mux.*x(:,1)*Xi1;
    Bx=Kxk./(Cx*Dx);
    Svx=Fz0*(beta(21)+beta(22)*dfz)*beta(6)*beta(2)*Xi1;
    Gxa=-1;
    Fx=(Dx.*sin(Cx*atan(Bx.*x(:,2)-Ex.*(Bx.*x(:,2)-atan(Bx.*x(:,2)))))+Svx)*Gxa;
end