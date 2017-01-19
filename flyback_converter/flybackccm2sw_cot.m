function [cost,rbp,rbn,rhn,rpn,rpp,rsn]=flybackccm2sw_cot(vinmin,vinmax,vomin,vomax,vin,vo,pi,eff,fs,ton,kdcm)
%constant on time flyback with two switch clamp in ccm/dcm(reverse current block) 
% usage:  flybackccm2sw_cot(37,57,3.3,5,48,5,6.5,0.89,2.5e5,1.2e-6,0.4);
%30% duty cycle seems about right
%rdson loss, dynamic swithing and metal routing are captured.
%vmaxp is maximum primary stress voltage. 
% lp is primary magnetizing inductance
% ton is fixed on time, fs is steady state ccm frequency
% kdcm is the ratio of io that any Iload less than that enters dcm for the given vin
% ideally to ensure ccm for the entire range of vin at kdcm*iload or lower, it should be solved for vinmax (toffmax_ccm)

close all;
cox=2.8e-9;  %f/mm2 assume the same for primary and secondary fets,2/3 is cgg and 1/3 is cgd
rwp85=40.88;    %Ohm*mm for 85V pldmos in magna chips
rwn85=15.85;     %for 85V nldmos in magna chips
rwn30=10.53;     %for 30V nldmos in global foundry
lp85=3.9e-3;      %drawn channel gate in mm
ln85=3e-3;
ln30=1.1e-3;
rmetal=0.04; %2squar metal routing for each drain and source (square shape device)
rmetals=0.02; %secondary FET metalization
ran85=0.9*0.18; %Ohm*mm2 for NLDMOS 85V, 0.9is a fudge factor from layout
rap85=1.3*0.58;   %PLDMOS 85V, 1.3 is a fudge factor from layout
ran30=0.06;     %NLDMOS 30V
cst85=0.03;     %cost in 85V process per mm2
cst30=0.05;     %cost in 30V process per mm2

vgate=5     %gate overdrive for dynamic switching
% each vo has a different reference design

po=eff*pi;
io=po/vo;

toffccm=1/fs-ton;  %in ccm
doff=toffccm*fs;
d=ton*fs;

n=ton*vin/(toffccm*vo);    
%in ccm and dcm, toffccm doesn't change by load,but increases with vin.
%however,  in dcm, fs does change with load (extended off time)
nmax=vinmin/vomax;  %limitted by primary shorted during off time
nmin=vinmax/(30-vomax); %limitted by sync fet stress voltage

toffmax_ccm=ton*vinmax/(n*vo);
dmin=ton/(ton+toffmax_ccm);
toffmin_ccm=ton*vinmin/(n*vo);
dmax=ton/(ton+toffmin_ccm);

%ensure CCM
% in general from gussian law: ipedccm*toffccm=io(ton+toffccm) (ipedccm=(ipk+ivl)/2), 
% and irips=n*vin/Lp*ton, 
% in ccm, transformer flux does not fully discharge every cycle, hence ivlsccm>0,

idcm=kdcm*io;   %boundry of entering dcm, any load greater than that is in ccm
ipeddcm=idcm/fs/toffccm;   %pedestal current at the boundry of entering dcm
% ipkdcms=2*ipeddcm (ivldcms=0)
irips=2*ipeddcm;   %delta I in secondary during off time in ccm/dcm
lp=vin*ton*n/irips;   %primary inductance, 
% for smaller vin it enters dcm at lower load current as irippls is smaller

ipedccm=io*(ton+toffccm)/toffccm;
ipksccm=ipedccm+0.5*irips;
ivlsccm=ipedccm-0.5*irips;   %in cot ivls does not come out negative

irmss=sqrt(doff/3*(ipksccm^2+ivlsccm^2+ipksccm*ivlsccm));  %seconday rms current
irmsp=sqrt(d/doff)*irmss/n;    %primary rms current

%in two switch clamp
vmaxp=vin;   %maximum primary fet voltage stress for a given vo
vmaxs=vin/n+vo;  %maximum secondary fet voltage stress for a given vo

ern85=0.5*rwn85*ln85*(2/3*cox*vgate^2+1/3*cox*(vmaxp+vgate)^2);   %switching-energy*Ohm for primary NFET, cgs=2/3cgg, cgd=1/3cgg
erp85=0.5*rwp85*lp85*(2/3*cox*vgate^2+1/3*cox*(vmaxp+vgate)^2);   %switching-energy*Ohm for primary PFET, cgs=2/3cgg, cgd=1/3cgg
ern30=0.5*rwn30*ln30*(2/3*cox*vgate^2+1/3*cox*(vmaxs+vgate)^2);   %switching-energy*Ohm for secondary NFET

r0=[2,2,1,0.5,0.03,1];    %initial rdson
r=[10,10,10,10,10,10];   
rlb=[0.001,0.001,0.001,0.001,0.001,0.001]; %lower bound rdson
rub=[10,10,10,10,10,10];    %upper bound rdson

% efficiency sweep
i=0;
for j=1:10
    i=i+1;
    effswp(i)=i*0.01+eff-0.07;
    [ri,costi]=fmincon(@(ri) mycost(ri,rap85,ran85,ran30,cst85,cst30),r0,[],[],[],[],rlb,rub, ...
    @(ri) mycon(ri,io,n,pi,effswp(i),ern30,ern85,erp85,fs,rmetal,rmetals,irmsp,irmss));
    costswp(i)=costi; costPi(i)=costi-ran30/ri(5)*cst30; costSi(i)=ran30/ri(5)*cst30;
    effi(i)=1-((rmetal+ri(1))*(io/n)^2+(rmetal+ri(2))*(io/n)^2+(rmetal+ri(3))*(io/n)^2+(2*rmetal+ri(4)+ri(6))*irmsp^2+(rmetals+ri(5))*irmss^2+ ...
     (ern30/ri(5)+ern85/ri(4)+erp85/ri(6))*fs)/pi;
    areaPi(i)=4*rap85/ri(1)+4*ran85/ri(2)+ran85/ri(3)+ran85/ri(4)+rap85/ri(6); areaSi(i)=ran30/ri(5);
    rbp_rpn(i)=ri(1)/ri(4); %normalized rdson, relatively constant across different efficiency points
    rbn_rpn(i)=ri(2)/ri(4);
    rhn_rpn(i)=ri(3)/ri(4);
    rpp_rpn(i)=ri(6)/ri(4);
    rsn_rpn(i)=ri(5)/ri(4);
    pbp_ppn(i)=(rmetal+ri(1))*(io/n)^2/((rmetal+ri(4))*irmsp^2+ern85/ri(4)*fs); %power loss distribution normalized vs total primary nfet loss
    pbn_ppn(i)=(rmetal+ri(2))*(io/n)^2/((rmetal+ri(4))*irmsp^2+ern85/ri(4)*fs);
    phn_ppn(i)=(rmetal+ri(3))*(io/n)^2/((rmetal+ri(4))*irmsp^2+ern85/ri(4)*fs);
    ppp_ppn(i)=((rmetal+ri(6))*irmsp^2+erp85/ri(6)*fs)/((rmetal+ri(4))*irmsp^2+ern85/ri(4)*fs);
    pdypp_ppn(i)=erp85/ri(6)*fs/((rmetal+ri(4))*irmsp^2+ern85/ri(4)*fs);
    pdypn_ppn(i)=ern85/ri(4)*fs/((rmetal+ri(4))*irmsp^2+ern85/ri(4)*fs);
end
figure('Name',sprintf('flyback-cot ccm efficiency sweep, vin=%d,vo=%0.1f,pi=%0.1f,d=%0.2f,n=%d,fs=%.2g',vin,vo,pi,d,n,fs),'NumberTitle','off');
hold on;
subplot(2,2,1);
plot(effi,costswp,'-b','linewidth',2);grid on;grid minor;
axis tight;title(sprintf('flyback-cot cost vs efficiency, vin=%d,vo=%0.1f,pi=%0.1f,d=%0.2f,n=%d,fs=%.2g',vin,vo,pi,d,n,fs));
ylabel('cost');xlabel('eff');
subplot(2,2,2);
plot(effi,areaPi,'-b',effi,areaSi,'-r','linewidth',2);grid minor;
axis auto;title('pri&sec area vs efficiency');
ylabel('mm2');xlabel('eff');
legend('primary','secondary',2);
subplot(2,2,3);
plot(effi,pbp_ppn,'-b',effi,pbn_ppn,'-r',effi,phn_ppn,'-k',effi,ppp_ppn,'-g',effi,pdypp_ppn,'-m',effi,pdypn_ppn,'-m','linewidth',2);grid minor;
axis auto;title('normalized total power per device vs efficiency');
ylabel('');xlabel('eff');
legend('Pbp/Ppn','Pbn/Ppn','Phn/Ppn','Ppp/Ppn','Pdy\_pp/Ppn','Pdy\_pn/Ppn',2);
subplot(2,2,4);
plot(effi,rbp_rpn,'-b',effi,rbn_rpn,'-r',effi,rhn_rpn,'-k',effi,rpp_rpn,'-g',effi,rsn_rpn,'-m','linewidth',2);grid minor;
axis auto;title('normalized Rsp vs efficiency');
ylabel('mm2');xlabel('eff');
legend('Rbp/Rpn','Rbn/Rpn','Rhn/Rpn','Rpp/Rpn','Rsn/Rpn',2);

[r,cost1]=fmincon(@(r) mycost(r,rap85,ran85,ran30,cst85,cst30),r0,[],[],[],[],rlb,rub,@(r) ...
    mycon(r,io,n,pi,eff,ern30,ern85,erp85,fs,rmetal,rmetals,irmsp,irmss));
rbp=r(1); wbp=rwp85/r(1);
rbn=r(2); wbn=rwn85/r(2);
rhn=r(3); whn=rwn85/r(3);
rpn=r(4); wpn=rwn85/r(4);
rpp=r(6); wpp=rwp85/r(6);
rsn=r(5); wsn=rwn30/r(5);

cost_total=cost1;
costP=cost1-ran30/r(5)*cst30;
costS=ran30/r(5)*cst30;
cstbg=(4*rap85/r(1)+4*ran85/r(2))*cst85;
csthp=(ran85/r(3))*cst85;
cstpfp=rap85/r(6)*cst85;
cstpfn=ran85/r(4)*cst85;

eff1=1-((rmetal+r(1))*(io/n)^2+(rmetal+r(2))*(io/n)^2+(rmetal+r(3))*(io/n)^2+(2*rmetal+r(4)+r(6))*irmsp^2+(rmetals+r(5))*irmss^2+ ...
    (ern30/r(5)+ern85/r(4)+erp85/r(6))*fs)/pi;
areaP=4*rap85/r(1)+4*ran85/r(2)+ran85/r(3)+ran85/r(4)+rap85/r(6);
areaS=ran30/r(5);
a_bp=rap85/r(1); a_bn=ran85/r(2); a_hn=ran85/r(3); a_pp=rap85/r(6); a_pn=ran85/r(4);

pwbp=(rmetal+r(1))*(io/n)^2; pwbn=(rmetal+r(2))*(io/n)^2; pwhp=(rmetal+r(3))*(io/n)^2; 
pwpp=(rmetal+r(6))*irmsp^2; pwpn=(rmetal+r(4))*irmsp^2; pwsn=(rmetals+r(5))*irmss^2;
pdpp=erp85/r(6)*fs; pdpn=ern85/r(4)*fs; pdsn=ern30/r(5)*fs;
cprip=rwp85*lp85*cox/r(6);cprin=rwn85*ln85*cox/r(4);csec=rwn30*ln85*cox/r(5);

f1=figure('Name','flyback summary','NumberTitle','off');
dat={'Rp_bg(Ohm)','Rn_bg(Ohm)','Rn_hp(Ohm)','Rp_pr(Ohm)','Rn_pr(Ohm)','Rn_sc(Ohm)','Rmetalp(Ohm)','Rmetalsec(Ohm)'; ...
    sprintf('%.2f',rbp), sprintf('%.2f',rbn), sprintf('%.2f',rhn), sprintf('%.2f',rpp), ...
    sprintf('%.2f',rpn), sprintf('%.3f',rsn), sprintf('%.3f',rmetal), sprintf('%.3f',rmetals); ...
    
    'Wp_bg(mm)','Wn_bg(mm)','Wn_hp(mm)','Wp_pr(mm)','Wn_pr(mm)','Wn_sc(mm)','Lp85(mm)','Ln85(mm)'; ...
    sprintf('%.2f',wbp), sprintf('%.2f',wbn), sprintf('%.2f',whn), sprintf('%.2f',wpp), ...
    sprintf('%.2f',wpn), sprintf('%.2f',wsn), sprintf('%.2g',lp85), sprintf('%.2g',ln85); ...
    
    'Ln30(mm)','A_bP/dev(mm2)','A_bN/dev(mm2)','A_hN(mm2)','A_priP(mm2)','A_priN(mm2)','A_pr(mm2)','A_sc(mm2)'; ...
    sprintf('%.2g',ln30), sprintf('%.2f',a_bp),sprintf('%.2f',a_bn),sprintf('%.2f',a_hn), ...
    sprintf('%.2f',a_pp),sprintf('%.2f',a_pn),sprintf('%.2f',areaP),sprintf('%.2f',areaS); ...
    
    'cost_total($)','cost_pr($)','cost_sc($)','cost_bridge','cost_hotplug','cost_prif_p','cost_prif_n','C_priP'; ...
    sprintf('%.2f',cost_total), sprintf('%.2f',costP), sprintf('%.2f',costS), sprintf('%.2f',cstbg),...
    sprintf('%.2f',csthp), sprintf('%.2f',cstpfp), sprintf('%.2f',cstpfn),sprintf('%.1g',cprip); ...

    'C_priN','C_sec','Pbridge_p','Pbridge_n','Photplug_n','Ppri_p','Ppri_n','Psec_n'; ...
    sprintf('%.1g',cprin), sprintf('%.1g',csec), sprintf('%.2f',pwbp), sprintf('%.2f',pwbn), ...
    sprintf('%.2f',pwhp), sprintf('%.2f',pwpp), sprintf('%.2f',pwpn), sprintf('%.2f',pwsn); ...
    
    'Pdyn_pri_p','Pdyn_pri_n','Pdyn_sec_n','io','irms_sec','Irip_sec','irms_pri','ipk_sec'; ...
    sprintf('%.2f',pdpp), sprintf('%.2f',pdpn), sprintf('%.2f',pdsn), sprintf('%.1f',io), ...
    sprintf('%.1f',irmss), sprintf('%.1f',irips), sprintf('%.1f',irmsp), sprintf('%.1f',ipksccm); ...
    
    'd','doff','dmax','dmin','n','fs','eff','vo'; ...
    sprintf('%.2f',d), sprintf('%.2f',doff), sprintf('%.2f',dmax), sprintf('%.2f',dmin), ...
    sprintf('%.1f',n),sprintf('%.2g',fs), sprintf('%.2f',eff1), sprintf('%.1f',vo); ...
    
    'pi','vin','vmax_pri','vmax_sec','Lpri','ton','toffccm','idcm'; ...
    sprintf('%.2f',pi),sprintf('%d',vin),sprintf('%.1f',vmaxp),sprintf('%.1f',vmaxs), ...
    sprintf('%.1e',lp),sprintf('%.2e',ton),sprintf('%.2e',toffccm),sprintf('%.2f',idcm)};
    
  
t=uitable('Parent',f1,'Data', dat,'ColumnName', {1,2,3,4,5,6,7,8},'RowName',{1,2,3,4,5,6,7,8}, ...
    'Units','Normalized','Position',[0.1 0.1 0.8 0.8],...
    'ColumnFormat',{'char','char','char'},'ColumnWidth','auto','FontUnits','normalized','FontSize',0.035);

%sweep vin for optimized flyback at full load (ccm)
i=0;
for j=1:21
    i=i+1;
    vinswp(i)=i+36;
    toffi(i)=ton*vinswp(i)/(n*vo);
    fsi(i)=1/(toffi(i)+ton);
    doffi(i)=toffi(i)*fsi(i);
    di(i)=ton*fsi(i); 

    ipedi(i)=io/fsi(i)/toffi(i);
    iripsi(i)=vinswp(i)*ton*n/lp;
    ipksi(i)=ipedi(i)+0.5*iripsi(i);
    ivlsi(i)=ipedi(i)-0.5*iripsi(i);   
    irmssi(i)=sqrt(doffi(i)/3*(ipksi(i)^2+ivlsi(i)^2+ipksi(i)*ivlsi(i)));
    irmspi(i)=sqrt(di(i)/doffi(i))*irmssi(i)/n;    
    effi(i)=1-((rmetal+r(1))*(io/n)^2+(rmetal+r(2))*(io/n)^2+(rmetal+r(3))*(io/n)^2+(2*rmetal+r(4)+r(6))*irmspi(i)^2+(rmetals+r(5))*irmssi(i)^2+ ...
    (ern30/r(5)+ern85/r(4)+erp85/r(6))*fsi(i))/pi;
end
figure('Name',sprintf('optimized ccm flyback-cot, n=%d,vo=%0.1f,po=%0.1f,ton=%.2g',n,vo,po,ton),'NumberTitle','off');
hold on;
subplot(1,3,1);
plot(vinswp,effi,'-b','linewidth',2);grid on;
axis auto;title(sprintf('optimized flyback-cot,eff vs vin,n=%d,vo=%0.1f,po=%0.1f,fs=%.2g',n,vo,po,fs));
ylabel('eff');xlabel('vin');
subplot(1,3,2);
plot(vinswp,di,'-r','linewidth',2);grid on;
axis auto;title('optimized flyback-cot,duty cycle vs vin');
ylabel('dutyC');xlabel('vin');
subplot(1,3,3);
plot(vinswp,fsi,'-r','linewidth',2);grid on;
axis auto;title('opt. flyback-cot,CCM freq vs vin');
ylabel('Fs\_ccm');xlabel('vin');


%sweep D for minimum cost in ccm flyback- constant on-time
j=0.02;
for i=1:20
    dk(i)=j+0.18;
    j=j+0.02;
    tonk(i)=dk(i)/fs;
    toffk(i)=1/fs-tonk(i);  %in ccm
    doffk(i)=toffk(i)*fs;
    nk(i)=tonk(i)*vin/(toffk(i)*vo);    
    toffmax_ccmk(i)=tonk(i)*vinmax/(nk(i)*vo);
    dmink(i)=tonk(i)/(tonk(i)+toffmax_ccmk(i));
    toffmin_ccmk(i)=tonk(i)*vinmin/(nk(i)*vo);
    dmaxk(i)=tonk(i)/(tonk(i)+toffmin_ccmk(i));
    ipeddcmk(i)=idcm/fs/toffk(i);
    iripsk(i)=2*ipeddcmk(i); 
    lpk(i)=vin*tonk(i)*nk(i)/iripsk(i);   
    ipedccmk(i)=io/fs/toffk(i);
    ipksccmk(i)=ipedccmk(i)+0.5*iripsk(i);
    ivlsccmk(i)=ipedccmk(i)-0.5*iripsk(i);   %in cot ivls does not come out negative
    irmssk(i)=sqrt(doffk(i)/3*(ipksccmk(i)^2+ivlsccmk(i)^2+ipksccmk(i)*ivlsccmk(i)));  
    irmspk(i)=sqrt(d/doffk(i))*irmssk(i)/nk(i);   
   
    [ri,costi]=fmincon(@(ri) mycost(ri,rap85,ran85,ran30,cst85,cst30),r0,[],[],[],[],rlb,rub, ...
    @(ri) mycon(ri,io,nk(i),pi,eff,ern30,ern85,erp85,fs,rmetal,rmetals,irmspk(i),irmssk(i)));
    costk(i)=costi; costPk(i)=costi-ran30/ri(5)*cst30; costSk(i)=ran30/ri(5)*cst30;
    effk(i)=1-((rmetal+ri(1))*(io/nk(i))^2+(rmetal+ri(2))*(io/nk(i))^2+(rmetal+ri(3))*(io/nk(i))^2+ ...
        (2*rmetal+ri(4)+ri(6))*irmspk(i)^2+(rmetals+ri(5))*irmssk(i)^2+(ern30/ri(5)+ern85/ri(4)+erp85/ri(6))*fs)/pi;
    
end
figure('Name',sprintf('flyback-cot ccm Duty cycle sweep vs cost\n vin=%d,vo=%0.1f,pi=%0.1f,eff=%0.2f,fs=%.2g',vin,vo,pi,eff,fs),'NumberTitle','off');
hold on;
subplot(2,2,1);
plot(dk,costk,'-r',dk,costPk,'-b',dk,costSk,'-g','linewidth',2);grid on;grid minor;
axis tight;title(sprintf('flyback-cot cost vs D\n vin=%d,vo=%0.1f,pi=%0.1f,eff=%0.2f,fs=%.2g',vin,vo,pi,eff,fs));
ylabel('cost');xlabel('duty cycle');
legend('cost','costP','costS','NorthEast');
subplot(2,2,2);
plot(dk,ipedccmk,'-b',dk,ipksccmk,'-r',dk,irmssk,'-m',dk,irmspk,'-g','linewidth',2);grid on;grid minor;
axis tight;title(sprintf('flyback-cot ac current vs D\n vin=%d,vo=%0.1f,pi=%0.1f,eff=%0.2f,fs=%.2g',vin,vo,pi,eff,fs));
ylabel('(A)');xlabel('duty cycle');
legend('ipedccmk','ipksccmk','irmssk','irmspk','North');
subplot(2,2,3);
plot(dk,nk,'-b','linewidth',2);grid on;grid minor;
line([dk(1),dk(20)],[nmin,nmin],'color','r','linewidth',2);
line([dk(1),dk(20)],[nmax,nmax],'color','r','linewidth',2);
axis tight;title(sprintf('flyback-cot number of turns vs D\n vin=%d,vo=%0.1f,pi=%0.1f,eff=%0.2f,fs=%.2g',vin,vo,pi,eff,fs));
ylabel('n');xlabel('duty cycle');
subplot(2,2,4);
plot(dk,lpk,'-b','linewidth',2);grid on;grid minor;
axis tight;title(sprintf('flyback-cot primary inductance vs D\n vin=%d,vo=%0.1f,pi=%0.1f,eff=%0.2f,fs=%.2g,idcm=%0.2f',vin,vo,pi,eff,fs,idcm));
ylabel('Lp');xlabel('duty cycle');


%sweep io in the optimized flyback- constant on-time, ccm/dcm
j=0;
for i=1:13
    ioe(i)=j+0.1;
    j=j+0.1;
    ipedccme(i)=ioe(i)*(ton+toffccm)/toffccm; %when ccm
    toffdcme(i)=ipeddcm*toffccm/ioe(i)-ton; %when dcm, toffdcme includes extended off time in dcm
    toffe(i)=max(toffccm,toffdcme(i));  %ccm/dcm combined
    fse(i)=1/(toffe(i)+ton);
    doffe(i)=toffccm*fse(i);
    de(i)=ton*fse(i);
    
    ipede(i)=max(ipedccme(i),ipeddcm);
    ipkse(i)=ipede(i)+0.5*irips;
    ivlse(i)=ipede(i)-0.5*irips;   %in cot ivls does not come out negative
    irmsse(i)=sqrt(doffe(i)/3*(ipkse(i)^2+ivlse(i)^2+ipkse(i)*ivlse(i)));  
    irmspe(i)=sqrt(de(i)/doffe(i))*irmsse(i)/n;   
   
    plosse(i)=(rmetal+r(1))*(ioe(i)/n)^2+(rmetal+r(2))*(ioe(i)/n)^2+(rmetal+r(3))*(ioe(i)/n)^2+ ...
        (2*rmetal+r(4)+r(6))*irmspe(i)^2+(rmetals+r(5))*irmsse(i)^2+(ern30/r(5)+ern85/r(4)+erp85/r(6))*fse(i);
    effe(i)=ioe(i)*vo/(ioe(i)*vo+plosse(i));

end

figure('Name',sprintf('Opt. flyback-cot load current sweep, vin=%d,vo=%0.1f,kdcm=%0.2f,iomax=%0.1f',vin,vo,kdcm,io),'NumberTitle','off');
hold on;
subplot(1,2,1);
plot(ioe,fse,'-r','linewidth',2);grid on;grid minor;
axis tight;title(sprintf('flyback-cot freq vs Io,\n vin=%d,vo=%0.1f,kdcm=%0.2f,iomax=%0.1f',vin,vo,kdcm,io));
ylabel('freq');xlabel('Io(A)');
subplot(1,2,2);
plot(ioe,effe,'-b','linewidth',2);grid on;grid minor;
axis tight;title(sprintf('flyback-cot eff (w/o iddq) vs Io,\n vin=%d,vo=%0.1f,kdcm=%0.2f,iomax=%0.1f',vin,vo,kdcm,io));
ylabel('%');xlabel('Io(A)');


end

function [c,ceq]=mycon(r,io,n,pi,eff,ern30,ern85,erp85,fs,rmetal,rmetals,irmsp,irmss)
c=(rmetal+r(1))*(io/n)^2+(rmetal+r(2))*(io/n)^2+(rmetal+r(3))*(io/n)^2+(2*rmetal+r(4)+r(6))*irmsp^2+(rmetals+r(5))*irmss^2+ ...
    (ern30/r(5)+ern85/r(4)+erp85/r(6))*fs-pi*(1-eff);
ceq=[];
end

function y=mycost(r,rap85,ran85,ran30,cst85,cst30)
y=(4*rap85/r(1)+4*ran85/r(2)+ran85/r(3)+ran85/r(4)+rap85/r(6))*cst85+ran30/r(5)*cst30;
end
