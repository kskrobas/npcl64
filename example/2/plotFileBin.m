fname='si111_7x7.laueb';

fin=fopen(fname,'r');

val=char(fread(fin,2,'uchar'));
val=char(fread(fin,1,'uchar'));


ksize=fread(fin,1,'uint64')

kmin=fread(fin,1,'double')
kstp=fread(fin,1,'double')
kmax=fread(fin,1,'double')

dataSize=fread(fin,1,'uint64')
nAtoms=fread(fin,1,'uint64')

# raf :  r - radiation  a - axis  f - format
raf=char(fread(fin,3,'uchar'))

lambda=fread(fin,1,'double')

## ignore 
fread(fin,1,'int');



val=char(fread(fin,4,'uchar'));


Kx=zeros(ksize,ksize,ksize);
Ky=zeros(ksize,ksize,ksize);

I=zeros(ksize,ksize,ksize);

##
for k=1:ksize
  for j=1:ksize
    for i=1:ksize
        
      I(i,j,k)=fread(fin,1,'double');      
      
    endfor
  endfor
  printf("\r   progress ....  %3.2f %%",k*100/ksize);
 endfor


fclose(fin);

K=kmin:kstp:kmax;



%  electron diffraction
%  lambda[nm] = sqrt(1.5/U)
%

 %[w,m]=max(max(max(I)));
 %surf(Kx(:,:,m),Ky(:,:,m),I(:,:,m));shading interp;view(0,-90)
 
 
 m=ceil(ksize*0.5);
 %m=11;
 
 %m
 ss=(m-2):(m+2);
 %I(ss,ss,m)=0;
 maxI=max(max(I(:,:,m)));
 I(:,:,m)/=maxI;
 I(ss,ss,m)/=1.25;
 
 figure(1)
 
 surf(K,K,I(:,:,m));shading interp;view(0,-90)

 axis square
 igray=1-gray;
 colormap(igray)
 
 
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
set(gca,'fontsize',12)


  hold on
  %r=1/lambda;
  r=0.5;
  alp=linspace(0,1,101)*2*pi;
  a=r*sin(alp);
  b=r*cos(alp);
  plot(a,b,'-r','linewidth',1)
  
  
   %img=flipud(imread('v00.png'));
   %image([0.25,1.25],[0.25,1.25],img)
 
  hold off
  
  if raf(1)=='X'
  title(['X-ray diffraction,  λ=' num2str(lambda)])
  else
  title(['Electron diffraction,  λ=' num2str(lambda)])
  end
  


%figure(2)
%
%Im=mean(I,3);
%
% surf(Kx(:,:,m),Ky(:,:,m),Im);shading interp;view(0,-90)
%
% axis square
% igray=1-gray;
% colormap(igray)
%
%xlabel('X','fontsize',14)
%ylabel('Y','fontsize',14)
%set(gca,'fontsize',12)
%


printf("\n");
