
unitCell_nOfAtoms=44;
unitCell_dim=[ 9.45614     0          0;
              -4.72807     8.18926    0;
               0           0       6.87830 ];
			
volUC=det(unitCell_dim);
aveDens=unitCell_nOfAtoms*1.1/volUC;

nOfSelAtoms=10;
volPart=2*3;
dr=0.125;				

datain=load('hap02.rdhl');
r=datain(:,1);
hss=volPart*sum(datain(:,2:end),2)/nOfSelAtoms;
gr=hss./(4*pi*r.^2.*dr)/aveDens;


plot([0 80],[1 1],'--k',r,gr,'-r','linewidth',2)
xlim([0 80])
set(gca,'fontsize',12)
xlabel('r','fontsize',14)
ylabel('g(r)','fontsize',14)
title('g(r)','fontsize',18)
grid
