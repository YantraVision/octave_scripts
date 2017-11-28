function [ APSKConst ] = APSK_Const_Gen( M, spc_eff )
%
%  function [ APSKConst ] = APSK_Const_Gen( M, spc_eff )
%  creates complex-valued symbols according to DVB-S2
%  it support 2 modulations 16 and 32 APSK.
%  Each circle radious depends on the spc_eff variable.
%
%     Inputs:
%       M: (scalar) modulation format (16 and 32 supported)
%       spc_eff: (scalar) variable to set the radious of the rings in the modulated
%       symbols.
%
%    Outputs:
%       APSKConst (vector): 16 or 32 APSK symbols (constellation)
%
%
% History:
%    Roberto Piazza   
%    University of Luxembourg
%    2014/06/13

switch(M)
    case 16
            
Circle(1)=4;
Circle(2)=12;
phase_shift(1)=0;
phase_shift(2)=0;
%/* Optimum capacity */
rho(1)=1;
if (spc_eff<0)
                rho(2)=2.57;
elseif(spc_eff<2.8)
                    rho(2)=3.15;
             
                elseif(spc_eff<3.1) 
                    rho(2)=2.85;
                
elseif(spc_eff<3.25)
                    rho(2)=2.75;
                
                elseif(spc_eff<3.4)
                    rho(2)=2.70;
                
elseif(spc_eff<3.57)    
                    rho(2)=2.60;
                else
                    rho(2)=2.57;
                end
        case 32
            Circle(1)=4;
Circle(2)=12;
Circle(3)=16;
phase_shift(1)      =0;
phase_shift(2)      =0;
phase_shift(3)      = pi/16; 
%/* DVB-S2 (Optimum capacity) */
            rho(1)      =1;
            if(spc_eff<0)           
                rho(2)=2.53;
rho(3)=4.30; 
            elseif(spc_eff<3.8) 
                    rho(2)=2.84;
rho(3)=5.27;
                elseif(spc_eff<4.1) 
                    rho(2)=2.72;
rho(3)=4.87;
                elseif(spc_eff<4.3) 
                    rho(2)=2.64;
rho(3)=4.64;
                elseif(spc_eff<4.45)
                    rho(2)=2.54;
rho(3)=4.33;
                else
                    rho(2)=2.53;
rho(3)=4.30;
            end
end



x=1;
for c_num=1:length(Circle)
   phase_offset = pi/Circle(c_num);
   for j = 1: Circle(c_num)
       phase = (2*j-1)*phase_offset + phase_shift(c_num);
       temp_APSKConst(x) = rho(c_num)*(cos(phase)+1i*sin(phase));
x=x+1;
   end
end
if (M==16)
APSKConst(1)=temp_APSKConst(6);
APSKConst(2)=temp_APSKConst(15);
APSKConst(3)=temp_APSKConst(9);
APSKConst(4)=temp_APSKConst(12);
APSKConst(5)=temp_APSKConst(5);
APSKConst(6)=temp_APSKConst(16);
APSKConst(7)=temp_APSKConst(10);
APSKConst(8)=temp_APSKConst(11);
APSKConst(9)=temp_APSKConst(7);
APSKConst(10)=temp_APSKConst(14);
APSKConst(11)=temp_APSKConst(8);
APSKConst(12)=temp_APSKConst(13);
APSKConst(13)=temp_APSKConst(1);
APSKConst(14)=temp_APSKConst(4);
APSKConst(15)=temp_APSKConst(2);
APSKConst(16)=temp_APSKConst(3);
elseif (M==32)
APSKConst(1)=temp_APSKConst(6);
APSKConst(2)=temp_APSKConst(7);
APSKConst(3)=temp_APSKConst(15);
APSKConst(4)=temp_APSKConst(14);
APSKConst(5)=temp_APSKConst(9);
APSKConst(6)=temp_APSKConst(8);
APSKConst(7)=temp_APSKConst(12);
APSKConst(8)=temp_APSKConst(13);
APSKConst(9)=temp_APSKConst(17);
APSKConst(10)=temp_APSKConst(19);
APSKConst(11)=temp_APSKConst(30);
APSKConst(12)=temp_APSKConst(28);
APSKConst(13)=temp_APSKConst(22);
APSKConst(14)=temp_APSKConst(20);
APSKConst(15)=temp_APSKConst(25);
APSKConst(16)=temp_APSKConst(27);
APSKConst(17)=temp_APSKConst(5);
APSKConst(18)=temp_APSKConst(1);
APSKConst(19)=temp_APSKConst(16);
APSKConst(20)=temp_APSKConst(4);
APSKConst(21)=temp_APSKConst(10);
APSKConst(22)=temp_APSKConst(2);
APSKConst(23)=temp_APSKConst(11);
APSKConst(24)=temp_APSKConst(3);
APSKConst(25)=temp_APSKConst(32);
APSKConst(26)=temp_APSKConst(18);
APSKConst(27)=temp_APSKConst(31);
APSKConst(28)=temp_APSKConst(29);
APSKConst(29)=temp_APSKConst(23);
APSKConst(30)=temp_APSKConst(21);
APSKConst(31)=temp_APSKConst(24);
APSKConst(32)=temp_APSKConst(26);
else
    APSKConst=temp_APSKConst;
end

APSKConst=APSKConst/sqrt(mean(abs(APSKConst).^2));
