%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to plot 1D-Spectra picture                            %
% Note:                                                                   %
%   You should run 'readStat_CH.m'/'readStat_HC.m' firstly                %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-12-17                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
function PlotSpectra
  global height PicYscale CBarXScale
  height=1;
  SpecDataDir='../../StatOut/';
  PicYscale=0.75;
  CBarXScale=0.6;
  
  VarStr='uu';DirStr='x';
  SpecFileDir=[SpecDataDir,'spec',DirStr,'_',VarStr,'.txt'];
  [SpecData,yplus,lamda,Retau]=GetSpecData(SpecFileDir);
  x=lamda/height*Retau; y=yplus;
  PlotSpec(x,y,SpecData,VarStr,DirStr,'inner',1);
  x=lamda/height; y=yplus/Retau;
  PlotSpec(x,y,SpecData,VarStr,DirStr,'outer',2);
end

function [SpecData,yplus,lamda,Retau]=GetSpecData(SpecFileDir)
  global height
  SpecData=load(SpecFileDir);
  nk=size(SpecData,1)-3;
  ny=size(SpecData,2)-1;
  Retau=SpecData(1,1);
  yplus=SpecData(1,2:end);
  waveNumber=SpecData(3:end-1,1);
  lamda=2*pi./waveNumber;

  SpecData=SpecData(3:end-1,2:end);
  for k=1:nk
    Ratio=waveNumber(k)*height;
    for m=1:ny
      SpecData(k,m)=SpecData(k,m)*Ratio;
    end
  end
end

function PlotSpec(x,y,SpecData,VarStr,DirStr,PicType,nFig)
  global PicYscale CBarXScale
  figure(nFig);hold on;
  [X,Y]=meshgrid(x,y);
  if(strcmp(VarStr,'uv')==0) 
    %pcolor(X,Y, SpecData');
    contourf(X,Y, SpecData',10,'k','LineWidth',0.1);
  else
    %pcolor(X,Y,-SpecData');  
    contourf(X,Y,-SpecData',10,'k','LineWidth',0.1);
  end
  shading interp;box on;
  if(strcmp(PicType,'inner')==1)
    set(gca,'XScale','log');set(gca,'YScale','log');
  else
    set(gca,'XScale','log');set(gca,'YScale','lin');  
  end
  set(gca,'FontSize',12,'FontName','Times New Roman');
  set(gca,'linewidth',0.7)
  set(gca,'TickDir','out','TickLength',[0.015,0.015]);

  set(gca,'XLim',[min(x),max(x)]);
  set(gca,'xMinorTick','on');
  if(strcmp(PicType,'inner')==1)
    xlabel(['\fontsize{16}\it\lambda_',DirStr,'^+']);
  else
    xlabel(['\fontsize{16}\it\lambda_',DirStr,'/h']);  
  end

  if(strcmp(PicType,'inner')==1)
    set(gca,'YLim',[1,max(y)]);
    ylabel('\fontsize{16}\ity^+')
    set(gca,'YMinorTick','on');
  else
    set(gca,'YLim',[0,1]);
    ylabel('\fontsize{16}\ity/h')
    set(gca,'YMinorTick','off');
  end

  ch=colorbar;
  set(ch,'FontSize',12,'FontName','Times New Roman');
  set(ch,'linewidth',0.6);
  if(strcmp(VarStr,'uv')==0)
    title(ch,['\fontsize{12}\itk_',DirStr,'E_{',VarStr,'}^+']);
  else
    title(ch,['\fontsize{12}-\itk_',DirStr,'E_{',VarStr,'}^+']);  
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
end