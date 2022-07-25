clc
close all
clear

% randf=rand(10^7,1); %%% this random vector is already generated and sa
% save randdata randf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randtype=0; %%%% random(non-reproducable results)
% randtype=1; %%%% reproducable results 

if randtype==1
    load randdata;
    count=1;
end

xmin=-10;  xmax=10;  dx=0.05;

%%% Interval
intv=xmin:dx:xmax;
N=max(size(intv));

%%%%%%%%%%%%%%%%%%%%%%%%%% objective function
xx=xmin:dx:xmax;
funx=fun(xx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cooling scheduale for VFSA  
nruns=1;
maxtemps=40;
nmov=3;
temp0=1;
decay=1.0;
t0=1;
nx=round(abs(xmax-xmin)./dx)+1;
w=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VFSA loop
temp=zeros(maxtemps,1);

for j=1:nruns
    ermod=0;
    if randtype==0
        xmod=sample1(xmin,dx,nx,xmin,xmax);
    elseif randtype==1
        [xmod,count]=sample1_reproducable(xmin,dx,nx,xmin,xmax,randf,count);
    else
        sprintf('randtype must be either 1 or zero')
        return
    end

    emod=fun(xmod);
    
    
    for jtemp=1:maxtemps
        temp(jtemp)=temp0.*exp(-decay.*(jtemp-1).^0.5);        
        tmp=t0.*exp(-decay.*(jtemp-1).^0.5);
        
                       
        for jmov=1:nmov
  
           if randtype==0 
               xtrial=walk(xmod,dx,xmin,xmax,tmp);
               
           else
               [xtrial,count]=walk_reproducable(xmod,dx,xmin,xmax,tmp,randf,count);
               
           end
            xtrial=round((xtrial-xmin)./dx).*dx+xmin;
            etrial=fun(xtrial);
            
            if etrial< emod
                emod=etrial;
                xmod=xtrial;                
            else
                arg=(etrial-emod)./temp(jtemp);
                if arg>1.e6
                    pde=0.001;
                else
                    pde=exp(-arg);
                end
                if randtype==0
                    if pde>rand
                        emod=etrial;
                        xmod=xtrial;
                    end
                else
                    if pde>randf(count)
                        emod=etrial;
                        xmod=xtrial;
                    end
                    count=count+1;
                end
            end              
            
        end %%%% end move
        ermod1(jtemp,j)=emod;
        m1(jtemp,j)=xmod;

    end
    qq=min(ermod1(:,j));
    best=find(qq==ermod1(:,j));
    if max(size(best))>1
        best=round(mean(best));
    end
    m2(j)=m1(best,j);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures
ss1=sum(sum(ermod1));
merror=mean(ermod1');

fs=16;

figure(1);
plot(xx,funx,'b','LineWidth',2);
title('Objective function','FontSize',fs);
xlabel('X','FontSize',fs);
ylabel('F(X)','FontSize',fs);
set(gca,'FontSize',fs)

figure(2);
plot(ermod1,'b','LineWidth',2);
hold on
xlim([1 maxtemps]);
xlabel('Iteration No.','FontSize',fs);
ylabel('Absolute error','FontSize',fs);
set(gca,'FontSize',fs)


figure(3);
plot(m1,'b','LineWidth',2);
hold on
xlim([1 maxtemps]);
xlabel('Iteration No.','FontSize',fs);
ylabel('X value','FontSize',fs);
set(gca,'FontSize',fs)













