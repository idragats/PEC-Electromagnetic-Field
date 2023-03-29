clc; clear;


%Different sigma, frequence = 1.0e9

 N=200;               %numper of cells 
 c=2.99792458e8;      %speed of light  
 fs=1.0e9;            %frequency
 T=1/fs;              %period
 lamda=c/fs;          %wavelength
 dx=lamda/20.0;       %space step
 dt=dx/c;             %time step 
 muz=4.0*pi*1.0e-7;   %permeability of free space
 epsz=1.0/(c*c*muz);  %permittivity of free space
 omega=2.0*pi*fs;
 omegadt=omega*dt;
 e=exp(1);
 nr=1/e; 

 Nt=30;                       %total number of periods to run simulation
 nmax=round(Nt*T/dt);         %total number of time steps


TIME=zeros(1,3);
delta=zeros(1,3);
SIGMA=zeros(1,3);

 epsr = ones(N,1);
 sigma = zeros(N,1);
 mur=ones(N,1);
 sim=zeros(N,1);
pl=1;
for profile=1:3
 
  for i=1:N
   
    if (profile==1)
      
       if (i>=100 && i<=160)
        epsr(i)=9;
        sigma(i)=0.01;
       end 
    end 
    
    if (profile==2)
      
       if (i>=100 && i<=160)
        epsr(i)=9;
        sigma(i)=0.1;
       end 
    end 
    
    if (profile==3)
      
       if (i>=100 && i<=160)
        epsr(i)=9;
        sigma(i)=1;
       end 
    end 
    
  end 
 
EE=zeros(N,nmax);
 

 ca=(1.0-(dt*sigma)./(2.0*epsz*epsr))./(1.0+(dt*sigma)./(2.0*epsz*epsr));
 cb=(dt/epsz./epsr/dx)./(1.0+(dt*sigma)./(2.0*epsz*epsr));

 da=(1.0-(dt*sim)./(2.0*muz*mur))./(1.0+(dt*sim)./(2.0*muz*mur));
 db=(dt/muz./mur/dx)./(1.0+(dt*sim)./(2.0*muz*mur));


 Ez=zeros(1,N);
 Hy=zeros(1,N);
 emaxz(1:N)=0;
 emaxy(1:N)=0;
 d=sqrt(10^6*2./(sigma.*mur.*omega));
 x=1:N;
 
 %TFSF boundary
 
 x1a=90*ones(1,N);
 x2a=90*ones(1,N);
 y1a=linspace(-3,3,N);
 y2a=linspace(-8e-3,8e-3,N);
 
 x1b=170*ones(1,N);
 x2b=170*ones(1,N);
 

 
 
j=1;


 for n=1:nmax
    
    Ez(1)=sin(omegadt*n); 
    
  
    for i=1:N-1
         Hy(i)=da(i)*Hy(i)-db(i)*(Ez(i+1)-Ez(i));
    end
    Hy(N)=0;
   
  
     for i=2:N
         Ez(i)=ca(i)*Ez(i)-cb(i)*(Hy(i)-Hy(i-1));
     end
     
    
     
     
     
     
     %mur's ABC first order
      Ez(1)=Ez(2); 
      Ez(N)=Ez(N-1);
    
    
     Ez(N)=0; %Pec
     
     
     emaxz=max(Ez,emaxz);
     emaxy=max(Hy,emaxy);
     
   EE(:,j)=emaxz;
     
  j=j+1;

     
 
  figure(pl)  
  subplot(2,1,1),plot(x,Ez,'r',x,emaxz,'g',x1a,y1a,'-black',x1b,y1a,'-black'),axis([0 N -3 3]);
  ylabel('Ez'); title(['n= ',num2str(round(n*dt*fs,1)),'ns  ','f= ',num2str(fs/10^9),'GHz  ','sigma= ',num2str(sigma(100))]);
  subplot(2,1,2),plot(x,Hy,'b',x,emaxy,'g',x2a,y2a,'-black',x2b,y2a,'-black'),axis([0 N -8e-3 8e-3]);
  ylabel('Hy'); xlabel('grid coordinate')
  pause(0.05)  
  
  
 end
pl= pl+1;  % for plots

time=zeros();
 k=1;
 for j=1:nmax
     for i=1:N
        if (EE(i,j)< nr)
            time(k)=j*dt*fs;
           k=k+1;
        end
     end
       
 end    
 
  TIME(profile)=time(1);
  SIGMA(profile)=sigma(100); 
end 

 
 
   for i =1:length(TIME)
         delta(i)=sqrt(2/(omega*TIME(i)*muz*SIGMA(i)));
        
   end
    
 figure(pl)   

 plot(SIGMA,delta,'-o') 
 hold on
 fplot(@(s) sqrt(2/(omega*muz*s)),[SIGMA(1) SIGMA(3)],'--')
 xlabel('sigma')
 ylabel('delta (micrometers)')
 title('skin depth')
 
 

%different frequenses sigma = 0.01

 
 N=200;                  %numper of cells 
 c=2.99792458e8;         %speed of light  
 FS=[1.0e9,3.0e9,6.0e9];   %frequences
 TIME=zeros(1,3);
 delta=zeros(1,3);
 pl=pl+1;
 for PROF=1:3
    
     if (PROF==1)
         fs=FS(1);
     end
      
     if (PROF==3)
         fs=FS(2);
     end
     
     if (PROF==9)
         fs=FS(3);
     end
 
 
 
 
 
 
 
 
 T=1/fs;              %period
 lamda=c/fs;          %wavelength
 dx=lamda/20.0;       %space step
 dt=dx/c;             %time step 
 muz=4.0*pi*1.0e-7;   %permeability of free space
 epsz=1.0/(c*c*muz);  %permittivity of free space
 omega=2.0*pi*fs;
 omegadt=omega*dt;
 e=exp(1);
 nr=1/e; 

 Nt=30;                       %total number of periods to run simulation
 nmax=round(Nt*T/dt);         %total number of time steps






 epsr = ones(N,1);
 sigma = zeros(N,1);
 mur=ones(N,1);
 sim=zeros(N,1);


 
  for i=1:N
   
   
      
       if (i>=100 && i<=160)
        epsr(i)=9;
        sigma(i)=0.01;
       end 
  
  end 
    
 
 
EE=zeros(N,nmax);
 

 ca=(1.0-(dt*sigma)./(2.0*epsz*epsr))./(1.0+(dt*sigma)./(2.0*epsz*epsr));
 cb=(dt/epsz./epsr/dx)./(1.0+(dt*sigma)./(2.0*epsz*epsr));

 da=(1.0-(dt*sim)./(2.0*muz*mur))./(1.0+(dt*sim)./(2.0*muz*mur));
 db=(dt/muz./mur/dx)./(1.0+(dt*sim)./(2.0*muz*mur));


 Ez=zeros(1,N);
 Hy=zeros(1,N);
 emaxz(1:N)=0;
 emaxy(1:N)=0;
 d=sqrt(10^6*2./(sigma.*mur.*omega));
 x=1:N;
 
 %TFSF boundary
 
 x1a=90*ones(1,N);
 x2a=90*ones(1,N);
 y1a=linspace(-3,3,N);
 y2a=linspace(-8e-3,8e-3,N);
 
 x1b=170*ones(1,N);
 x2b=170*ones(1,N);
 

 
 
j=1;


 for n=1:nmax
    
    Ez(1)=sin(omegadt*n); 
    
  
    for i=1:N-1
         Hy(i)=da(i)*Hy(i)-db(i)*(Ez(i+1)-Ez(i));
    end
    Hy(N)=0;
   
  
     for i=2:N
         Ez(i)=ca(i)*Ez(i)-cb(i)*(Hy(i)-Hy(i-1));
     end
     
    
     
     
     
     
     %mur's ABC first order
      Ez(1)=Ez(2); 
      Ez(N)=Ez(N-1);
    
    
     Ez(N)=0; %Pec
     
     
     emaxz=max(Ez,emaxz);
     emaxy=max(Hy,emaxy);
     
   EE(:,j)=emaxz;
     
  j=j+1;

     
 
  figure(pl)  
  subplot(2,1,1),plot(x,Ez,'r',x,emaxz,'g',x1a,y1a,'-black',x1b,y1a,'-black'),axis([0 N -3 3]);
  ylabel('Ez'); title(['n= ',num2str(round(n*dt*fs,1)),'ns  ','f= ',num2str(fs/10^9),'GHz  ','sigma= ',num2str(sigma(100))]);
  subplot(2,1,2),plot(x,Hy,'b',x,emaxy,'g',x2a,y2a,'-black',x2b,y2a,'-black'),axis([0 N -8e-3 8e-3]);
  ylabel('Hy'); xlabel('grid coordinate')
  pause(0.05)  
  
  
 end
 

 
  
 
 
time=zeros();
 k=1;
 for j=1:nmax
     for i=1:N
        if (EE(i,j)< nr)
            time(k)=j*dt;
           k=k+1;
           
        end
     end
       
 end    
 
 pl=pl+1;
   TIME(PROF)=time(1)*FS(PROF);
 
 
 end
 
 
   for i =1:length(TIME)
         delta(i)=sqrt(1/(pi*FS(i)*TIME(i)*muz*sigma(100)));
   end
    
   
 figure(pl)   
 
 plot(FS,delta,'o-') 
 hold on 
 fplot(@(f) sqrt(1/(pi*f*muz*0.01)),[FS(1) FS(3)],'--')
 xlabel('Frequence Hz')
 ylabel('delta (micrometers)')
 title('skin depth')

 