%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to read and calculate 2D-Spectra data                 %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-12-17                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&% 
clc;clear;

nxc=1024;
nzc=512;
nySpec2D=9;
NEnergySpec2D=5;
real_prec='real*8';

iTimeSet=0;
VarStr={'uu','vv','ww','pp','uv','cc'};

dir_statIn = '../../CFD/Results/';
dir_statOut='../../StatOut/';

%% ========== Normally no need to change anything below ==========
if(strcmp(real_prec,'real*4')==1) 
  real_byte=4;
elseif(strcmp(real_prec,'real*8')==1)
  real_byte=8;
else
  error('readSpec2D: real_prec wrong');
end
if ( exist(dir_statOut,'dir') ==false )
  mkdir(dir_statOut);
  fprintf( '%s\n\n',['crete directory: ',dir_statOut,'  sucessfully'] );
end
nxh=nxc/2; nxhp=nxh+1;
nzh=nzc/2; nzhp=nzh+1;
dir_output=dir(fullfile(dir_statIn,'Spec2D*'));
file_names={dir_output.name};
file_num=length(file_names);
if(file_num>0)
  for k=1:file_num
    datapath = cell2mat(file_names(k));
    if(str2double(datapath(7:16))<=iTimeSet);continue;end;
    datapath = [dir_statIn,cell2mat(file_names(k))];
    fid=fopen(datapath,'r');
    fseek(fid, 0, 'eof'); 
    totalbyte1=NEnergySpec2D*nySpec2D*nxc*nzc*real_byte;
    totalbyte2=ftell(fid);
    if(totalbyte1~=totalbyte2) 
      error('readSpec2D: File Byte Wrong');
    end
    fclose(fid);      
  end
  for NE=1:NEnergySpec2D
    file_ave=0;
    Spectra2DTemp=zeros(nxc,nySpec2D,nzc);
    offset=(NE-1)*nySpec2D*nxc*nzc*real_byte;
    for k=1:file_num
      datapath = cell2mat(file_names(k));
      if(str2double(datapath(7:16))<=iTimeSet);continue;end;
      datapath = [dir_statIn,cell2mat(file_names(k))];
      file_ave=file_ave+1;
      fid=fopen(datapath,'r');
      fseek(fid, offset, 'bof');
      SpecData=fread(fid,nySpec2D*nxc*nzc,real_prec);
      Spectra2DTemp=Spectra2DTemp+reshape(SpecData,size(Spectra2DTemp));
      fclose(fid);
      disp( [cell2mat(VarStr(NE)),' read:   ',datapath,'  sucessfully'] );
    end
    Spectra2DTemp=Spectra2DTemp/file_ave;
    Spectra2D=zeros(nxhp,nzhp,nySpec2D);
    for jk=1:nySpec2D
      Spectra2D(1,   1,   jk)= Spectra2DTemp(1,   jk,1   );
      Spectra2D(1,   nzhp,jk)= Spectra2DTemp(1,   jk,nzhp);
      Spectra2D(nxhp,1,   jk)= Spectra2DTemp(nxhp,jk,1   );
      Spectra2D(nxhp,nzhp,jk)= Spectra2DTemp(nxhp,jk,nzhp);
      m=1;
      for k=1:nzh-1
        k1=k+1; k2=nzc+2-k1;
        Spectra2D(m,k1,jk)=2.0*(Spectra2DTemp(m,jk,k1)+Spectra2DTemp(m,jk,k2));
      end  
      k=1;
      for m=1:nxh-1
        m1=m+1; m2=nxc+2-m1;
        Spectra2D(m1,k,jk)=2.0*(Spectra2DTemp(m1,jk,k)+Spectra2DTemp(m2,jk,k));
      end
      for k=1:nzh-1
        k1=k+1; k2=nzc+2-k1;
        for m=1:nxh-1
          m1=m+1; m2=nxc+2-m1;
          Spectra2D(m1,k1,jk)=4.0*(Spectra2DTemp(m1,jk,k1)+Spectra2DTemp(m2,jk,k1)+Spectra2DTemp(m1,jk,k2)+Spectra2DTemp(m2,jk,k2));
        end
      end  
    end
    writename=[dir_statOut,'Spec2D_',cell2mat(VarStr(NE))];
    fid=fopen(writename,'w');
    fseek(fid,0,'bof');
    fwrite(fid,Spectra2D,real_prec);
    fclose(fid);
  end
end