%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to plot 2D-Spectra picture                            %
% Note:                                                                   %
%   You need to run readSpec2D.m firstly.                                 %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-12-17                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&% 
clc;clear;

xlx=8*pi;
zlz=3*pi;
nxc=1024;
nzc=512;
nySpec2D=9;
real_prec='real*8';
utau=0.0428571;
height=1;
FileStr='../../StatOut/Spec2D_uu';
nyIndex=5;
PicYscale=0.75;
CBarXScale=0.6;
  
%% ========== Normally no need to change anything below ==========
nxh=nxc/2; nxhp=nxh+1;
nzh=nzc/2; nzhp=nzh+1;
if(strcmp(real_prec,'real*4')==1) 
  real_byte=4;
elseif(strcmp(real_prec,'real*8')==1)
  real_byte=8;
else
  error('readSpec2D: real_prec wrong');
end
fid=fopen(FileStr,'r');
fseek(fid, 0, 'eof');
totalbyte1=nxhp*nzhp*nySpec2D*real_byte;
totalbyte2=ftell(fid);
if(totalbyte1~=totalbyte2) 
  error('readSpec2D: File Byte Wrong');
end
fseek(fid,0,'bof');
SpecData=fread(fid,nySpec2D*nxhp*nzhp,real_prec');
fclose(fid);

SpecData=reshape(SpecData,[nxhp,nzhp,nySpec2D]);
SpecData=SpecData(:,:,nyIndex)/utau^2;
SpecData=SpecData(2:end-1,2:end-1);

waveNumberX=(1:nxh-1)*2*pi/xlx;
lamdaX=2*pi./waveNumberX';
waveNumberZ=(1:nzh-1)*2*pi/zlz;
lamdaZ=2*pi./waveNumberZ';
for m=1:nzh-1
  for k=1:nxh-1
    SpecData(k,m)=SpecData(k,m)*waveNumberX(k)*waveNumberZ(m);  
  end
end
VarStr=FileStr(end-1:end);
if(strcmp(VarStr,'uv')==1)
  SpecData=-SpecData; 
end
[X,Y]=meshgrid(lamdaX,lamdaZ);

%pcolor(X,Y, SpecData');
contourf(X,Y, SpecData',10,'k','LineWidth',0.1);
shading interp;box on;
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'FontSize',12,'FontName','Times New Roman');
set(gca,'linewidth',0.7)
set(gca,'TickDir','out','TickLength',[0.015,0.015]);
set(gca,'xMinorTick','on');
set(gca,'yMinorTick','on');
xlabel(['\fontsize{16}\it\lambda_x/h']);
ylabel(['\fontsize{16}\it\lambda_z/h']);

ch=colorbar;
set(ch,'FontSize',12,'FontName','Times New Roman');
set(ch,'linewidth',0.6);
if(strcmp(VarStr,'uv')==1)
  title(ch,['\fontsize{12}\itk_xk_zE_{',VarStr,'}^+']);
else
  title(ch,['\fontsize{12}\it-k_xk_zE_{',VarStr,'}^+']); 
end
set(ch,'TickDir','in','TickLength',0.02);
color=[255 255 255;185 255 255;155 255 255;
        45 255 255; 35 255 220;135 255 120;
       236 255  19;255 201   0;255 134   0;
       255  67   0;255  0    0]/255;
colormap(color); %get current colormap
if(strcmp(VarStr,'uv')==1)
 set(ch,'YLim',[0,-min(min(SpecData))]);
end
axpos=get(gca,'Position');
axpos(4)=PicYscale*axpos(4);
ch.Position(3)=CBarXScale*ch.Position(3);
ch.Position(4)=axpos(4);
set(gca,'Position',axpos);
hold off;