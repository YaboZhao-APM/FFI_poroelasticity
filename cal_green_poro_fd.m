function [AMARr,G,U_obs,bdata]=cal_green_poro_fd(fault_xyz,insardata1,p,a,k)
    distance=[];theta=[];AMAR=[];
    AMARr=[];
%     insardata=insardata3;
    obs_loc_xy=insardata1(:,1:2);
    obs_ve=insardata1(:,4);
    obs_vn=insardata1(:,5);
    obs_vz=insardata1(:,6);

  %
data_local=fault_xyz(:,1:3);

    %%%=================================
    for i=1:size(data_local,1)
        distancea=sqrt((data_local(i,1)-obs_loc_xy(:,1)).^2+(data_local(i,2)-obs_loc_xy(:,2)).^2)./1000;
        a1=obs_loc_xy(:,2)-data_local(i,2);a2=obs_loc_xy(:,1)-(data_local(i,1));
        thetaa=atan2d(a1,a2);
        distance=[distance distancea];
        theta=[theta thetaa];
    end
    %make mar A
    [x,y]=size(distance);
    for j=1:y
        gffile=sprintf('grn_%d',ceil(fault_xyz(j,3)));
        floder=sprintf('green_%3.2f_%3.2f_%E\\',p,a,k);
        pathname='C:\Users\lenovo\Desktop\model_poro\';
        GF=load([pathname floder gffile]);
%                 GF=load(gffile);
            delta=fault_xyz(j,5);%璁惧???捐?
         cD2=cosd(2*delta); % cos(2*dip)
     sD2=sind(2*delta); % sin(2*dip)
     cD=cosd(delta);    % cos(dip)
     sD=sind(delta);    % sin(dip)
%         gf12=GF(:,[1,2:4]);gf32=GF(:,[1,5:7]);gf22=GF(:,[1,8:10]);gf33=GF(:,[1,11:13]);
for i=1:x
         [sGFv,fGFv,wpGFv,tGFv]=InterpolationGF(GF,distance(i,j));
         strike=fault_xyz(j,4);
            theta_s=theta(i,j)+(90-strike);
             cS2=cosd(2*theta_s); % cos(2*theta)
             sS2=sind(2*theta_s); % sin(2*theta)           
             cS=cosd(theta_s);    % cos(theta)
             sS=sind(theta_s);    % sin(theta)
            traz1=sGFv(1,1)*sS2*sD-sGFv(1,2)*cS*cD;
            traz2=0.5*(sGFv(1,4)-sGFv(1,3)+sGFv(1,1)*cS2)*sD2 + sGFv(1,2)*sS*cD2;
            trar1=sGFv(2,1)*sS2*sD-sGFv(2,2)*cS*cD;
            tr1= tGFv(2,1)*sS2*sD-tGFv(2,2)*cS*cD;
            tr2=0.5*tGFv(2,1)*cS2*sD2 + tGFv(2,2)*sS*cD2;
            trar2=0.5*(sGFv(2,4)-sGFv(2,3)+sGFv(2,1)*cS2)*sD2 + sGFv(2,2)*sS*cD2;
            trathe1=sGFv(3,1)*cS2*sD+sGFv(3,2)*sS*cD;
            tt1=tGFv(3,1)*cS2*sD+tGFv(3,2)*sS*cD;
            trathe2=-0.5*sGFv(3,1)*sS2*sD2 + sGFv(3,2)*cS*cD2;
            tt2=-0.5*tGFv(3,1)*sS2*sD2 + tGFv(3,2)*cS*cD2;
            trar11=tr1+trar1;
            trar22=tr2+trar2;
            trathe11=trathe1+ tt1;
             trathe22=trathe2+tt2;
            %zrt2esup
            theta_ss=theta(i,j);
            COS=cosd(theta_ss);
            SIN=sind(theta_ss);
            A=[trar11*COS-trathe11*SIN  trar22*COS-trathe22*SIN fGFv(2,1)*COS;...
                trar11*SIN+trathe11*COS trar22*SIN+trathe22*COS fGFv(2,1)*SIN;...
                -1*traz1 -1*traz2 -fGFv(1,1)].*fault_xyz(j,6).*fault_xyz(j,7);
            xyz2los1=[obs_ve(i) -1*obs_vn(i) obs_vz(i)];
            A2=xyz2los1*A;
            AMAR=[AMAR;A2];
        end
        AMARr=[AMARr AMAR];
        AMAR=[];
    end
%     for i=1:length(obs_ve)
%         xyz2los1=[obs_ve(i); -1*obs_vn(i) ;obs_vz(i)];
%         xyz2los=[xyz2los ;xyz2los1];
%     end
%     AMARr=xyz2los.*AMARr;%xyz2los
%     AMARr=reshape(sum(reshape(AMARr,3,[])),[],size(AMARr,2)) ;
U_obs=insardata1(:,3);
 h2 = size(U_obs,1);
 rms_insar = ones(h2,1);      % uniform weighting for InSAR/AZO data
w_insar = calc_weight_insar_error(rms_insar);
bdata = w_insar * U_obs;  %cm
G = w_insar * AMARr; %cm
G1=G(:,1:3:end);AMARr1=AMARr(:,1:3:end);
G2=G(:,2:3:end);AMARr2=AMARr(:,2:3:end);
G3=G(:,3:3:end);AMARr3=AMARr(:,3:3:end);
G=[G1 G2 G3];AMARr=[AMARr1 AMARr2 AMARr3];
    
    
end



%%%end
    function [sGFv,fGFv,wpGFv,tGFv]=InterpolationGF(GF,Dis) %2D =(3,4)=(z,r,theta;12,32,22,33)
   [n,m]=size(GF);
   for j=2:m
       for i=1:n-1
           if Dis>=GF(i,1) && Dis<GF(i+1,1)
               break;
           end
       end
       iGF(j)=( GF(i+1,j)-GF(i,j) )/( GF(i+1,1)-GF(i,1) )*(Dis-GF(i,1))+ GF(i,j);
       x= GF(i+1,j)-GF(i,j);
       y=GF(i+1,1)-GF(i,1);
       z=Dis-GF(i,1);
       t=GF(i,j);
   end
   sGFv(:,1)=iGF(2:4);
   sGFv(:,2)=iGF(5:7);
   sGFv(:,3)=iGF(8:10);
   sGFv(:,4)=iGF(11:13);
   fGFv(:,1)=iGF(14:16);
   wpGFv(:,1)=iGF(17:19);
   tGFv(:,1)=iGF(20:22);
   tGFv(:,2)=iGF(23:25);
   

end