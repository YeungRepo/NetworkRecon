%clear all 
load x4_ey_ucsb_dat.mat

close all
 
offset = 0; 
mint_1 = 1;
maxt_1 = 49;
mint_2 = 49+offset+1;
maxt_2 = 150;
ara_first_1uM = [7,8,9];
ara_first_1mM = [10,11,12];
hsl_1nM_first = [19,20,21];
hsl_1uM_first = [22,23,24];
ara_then_hsl_1uM = [7,8,9];
ara_then_hsl_1mM = [10,11,12];
hsl_1nM_then_ara = [19,20,21];
hsl_1uM_then_ara = [22,23,24];
nc_1 = [34,35,36];
nc_2 = [34,35,36];
evd_1 = [28,29,30];
evd_2 = [28,29,30];
t = 8:10:(maxt_2-offset)*10;
t = t/60';





%offset = OD600_2(3,:) - OD600_2(1,:);
%od_bg_2 = repmat(offset,size(OD600_2,1),1);
%od = [OD600;max(OD600_2-odbg_2,0)];
od = [table2array(OD600);table2array(OD600_2)];
cfp = [table2array(CFPGain100_1);table2array(CFPGain100_2)];
yfp = [table2array(YFPGain61_1);table2array(YFPGain61_2)];
rfp = [table2array(RFPGain61_1);table2array(RFPGain61_2)];

dilfac = 20/500;
sm = 1/dilfac

for i = 1:size(cfp,1)
    for j = 1:size(cfp,2)
        cfpn(i,j) = cfp(i,j)/od(i,j);
        yfpn(i,j) = yfp(i,j)/od(i,j);
        rfpn(i,j) = rfp(i,j)/od(i,j);
        
    end
end


%rfp = [max(RFPGain150-rfpbg_1,0);max(RFPGain150_2-rfpbg_2,0)];
%yfpbg = max(max(yfp(mint_1:maxt_1,[nc_2]),[],1),[],2);
%cfpbg = max(max(cfp(mint_1:maxt_1,[nc_2]),[],1),[],2);

%yfpbg = repmat(mean(yfp(:,[nc_2]),2),1,size(yfp,2)); 

cfpbg_1 = repmat(max(cfpn(mint_1:maxt_1,evd_1),[],2),1,size(CFPGain61_1,2));
cfpbg_2 = repmat(max(cfpn(maxt_1+1:end,evd_1),[],2),1,size(CFPGain61_1,2));
cfpbg = [cfpbg_1;cfpbg_2]%repmat(mean(cfp(:,[nc_2]),2),1,size(cfp,2)); 


yfpbg_1 = repmat(max(yfpn(mint_1:maxt_1,evd_1),[],2),1,size(YFPGain61_1,2));
yfpbg_2 = repmat(max(yfpn(maxt_1+1:end,evd_1),[],2),1,size(YFPGain61_1,2));
yfpbg = [yfpbg_1;yfpbg_2];%repmat(mean(cfp(:,[nc_2]),2),1,size(cfp,2)); 

rfpbg_1 = repmat(mean(rfpn(mint_1:maxt_1,evd_1),2),1,size(RFPGain61_1,2));
rfpbg_2 = repmat(mean(rfpn(maxt_1+1:end,evd_1),2),1,size(RFPGain61_1,2));
rfpbg = [rfpbg_1;rfpbg_2];%repmat(mean(cfp(:,[nc_2]),2),1,size(cfp,2)); 



%Turned off background subtraction 
cfpn = max(cfpn,0);
yfpn = max(yfpn,0);
rfpn = max(rfpn,0);

%cfpn = max(cfpn-cfpbg,0);
%yfpn = max(yfpn-yfpbg,0);
%rfpn = max(rfpn-rfpbg,0);


%yfp = max(yfp-yfpbg,0);
%rfp = max(rfp-rfpbg,0);
%cfp = max(cfp-cfpbg,0);

odbg = 0;% min(min(od(mint_1:maxt_2,nc_2),[],1),[],2);



%od = max(od-odbg,0);


climit = 5e4;
ylimit = 2e3;
rlimit = 8e4;
%cfpbg = repmat(max([cfpn(mint_1:maxt_1,[nc_1]); cfpn(mint_2:maxt_2,[nc_2])],[],2),1,size(cfpn,2));


ct = cfpn;%max(cfpn(1:maxt_2,:)-0,0);
yt = yfpn;%max(yfpn(1:maxt_2,:)-0,0);
rt = rfpn;%max(rfpn(1:maxt_2,:)-0,0);



figc = figure(1);
figy = figure(2);
figr = figure(3);


% - - - - - - - plot yfp  - - - - - 
dilfac = 20/500;
sm = 1/dilfac
channel = yt;
colovec = [.5 .5 0];
catind = {ara_first_1uM,ara_then_hsl_1uM};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

xlimits =[t(mint_1) t(maxt_2-offset)];
figy = make_pretty_plot(t,y_plot,figy,colovec,[0 ylimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white')

catind = {hsl_1nM_first,hsl_1nM_then_ara};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

figy = make_pretty_plot(t,y_plot,figy,colovec.^2,[0 ylimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white','location','northwest','FontSize',20)


% - - - plot cyan 
channel = ct;
colovec = [0 0 0.9];
catind = {ara_first_1uM,ara_then_hsl_1uM};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

%xlimits =[t(mint_1) t(maxt_2)];
figc = make_pretty_plot(t,y_plot,figc,colovec,[0 climit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white')

catind = {hsl_1nM_first,hsl_1nM_then_ara};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

figc = make_pretty_plot(t,y_plot,figc,colovec.^16,[0 climit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white','location','northwest','FontSize',20)

%- - - - plot red 

channel = rt;
colovec = [.9 0 0 ];
catind = {ara_first_1uM,ara_then_hsl_1uM};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

%xlimits =[t(mint_1) t(maxt_2)];
figr = make_pretty_plot(t,y_plot,figr,colovec,[0 rlimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white')

catind = {hsl_1nM_first,hsl_1nM_then_ara};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

figr = make_pretty_plot(t,y_plot,figr,colovec.^16,[0 rlimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white','location','northwest','FontSize',20)




% - - - arabinose at 1 mM , HSL 1 uM -  - - - - 


figc = figure(4);
figy = figure(5);
figr = figure(6);




% - - - - - - - plot yfp  - - - - - 
dilfac = 20/500;
sm = 1/dilfac
channel = yt;
colovec = [.5 .5 0];
catind = {ara_first_1mM,ara_then_hsl_1mM};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

%xlimits =[t(mint_1) t(maxt_2)];
figy = make_pretty_plot(t,y_plot,figy,colovec,[0 ylimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white')

catind = {hsl_1uM_first,hsl_1uM_then_ara};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

figy = make_pretty_plot(t,y_plot,figy,colovec.^2,[0 ylimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white','location','northwest','FontSize',20)


% - - - plot cyan 
channel = ct;
colovec = [0 0 0.9];
catind = {ara_first_1mM,ara_then_hsl_1mM};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

%xlimits =[t(mint_1) t(maxt_2)];
figc = make_pretty_plot(t,y_plot,figc,colovec,[0 climit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white')

catind = {hsl_1uM_first,hsl_1uM_then_ara};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

figc = make_pretty_plot(t,y_plot,figc,colovec.^16,[0 climit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white','location','northwest','FontSize',20)

%- - - - plot red 

channel = rt;
colovec = [.9 0 0 ];
catind = {ara_first_1mM,ara_then_hsl_1mM};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

%xlimits =[t(mint_1) t(maxt_2)];
figr = make_pretty_plot(t,y_plot,figr,colovec,[0 rlimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white')

catind = {hsl_1uM_first,hsl_1uM_then_ara};
y_plot = [channel(mint_1:maxt_1,catind{1});channel(mint_2:maxt_2,catind{2})];

figr = make_pretty_plot(t,y_plot,figr,colovec.^16,[0 rlimit],xlimits);
h_leg = legend('','')
set(h_leg,'EdgeColor','white','location','northwest','FontSize',20)




%Y1 = ct(:,[ara_first_1nM,ara_first_1uM,hsl_1nM_first,hsl_1uM_first]); 
Y1 = yt(1:maxt_2,[ara_first_1uM,ara_first_1mM,hsl_1nM_first,hsl_1uM_first]); 
Y1 = Y1/max(max(Y1));

Y2 = rt(1:maxt_2,[ara_first_1uM,ara_first_1mM,hsl_1nM_first,hsl_1uM_first]); 
Y2 = Y2/max(max(Y2));

for col_index = 1:size(Y1,2) 
    Y1(:,col_index) = smooth(Y1(:,col_index),50); 
    Y2(:,col_index) = smooth(Y2(:,col_index),50); 
end
for row_index = 1:size(Y1,1)
    for col_index = 1:size(Y1,2)
        Y1(row_index,col_index) = max(0.0 ,Y1(row_index,col_index)-0.05);
        Y2(row_index,col_index) = max(0.0 ,Y2(row_index,col_index)-0.05);
    end
end


nM_input_train = 1.0/3.0*ones(maxt_1,1);
uM_input_train = ones(maxt_1,1);

nM_input_train_2 = 1.0/3.0*ones(maxt_2-maxt_1,1);
uM_input_train_2 = ones(maxt_2-maxt_1,1);

zeros_input_train =zeros(maxt_1,1);
U1 = [[repmat(nM_input_train,1,3) ; repmat(nM_input_train_2,1,3)] , [repmat(uM_input_train,1,3); repmat(uM_input_train_2,1,3)],  [repmat(zeros_input_train,1,3); repmat(nM_input_train_2,1,3)] ,[repmat(zeros_input_train,1,3); repmat(uM_input_train_2,1,3)] ];
U2 = [[repmat(zeros_input_train,1,3);repmat(nM_input_train_2,1,3)] ,[repmat(zeros_input_train,1,3);repmat(uM_input_train_2,1,3)], [repmat(nM_input_train,1,3) ; repmat(nM_input_train_2,1,3)] ,[repmat(uM_input_train,1,3); repmat(uM_input_train_2,1,3)]  ];

%Prepend_Zeros_Window = 0; 
%Y1 = [zeros(Prepend_Zeros_Window,size(Y1,2));Y1];
%Y2 = [zeros(Prepend_Zeros_Window,size(Y1,2));Y2];
%U1 = [zeros(Prepend_Zeros_Window,size(Y1,2));U1];
%U2 = [zeros(Prepend_Zeros_Window,size(Y1,2));U2];

c1r1 = 1;
c1r2 = 2; 
c1r3 = 3; 
c2r1 = 4;
c2r2 = 5;
c2r3 = 6;

c3r1 = 7; 
c3r2 = 8;
c3r3 = 9;

c4r1 = 10;
c4r2 = 11; 
c4r3 = 12; 

d1 = iddata([Y1(:,c1r1), Y2(:,c1r1)],[Y1(:,c1r1), Y2(:,c1r1), U1(:,c1r1),U2(:,c1r1)],'Ts',1,'TimeUnit','hours');
d2 = iddata([Y1(:,c1r2), Y2(:,c1r2)],[Y1(:,c1r2), Y2(:,c1r2), U1(:,c1r2),U2(:,c1r2)],'Ts',1,'TimeUnit','hours');
d3 = iddata([Y1(:,c1r3), Y2(:,c1r3)],[Y1(:,c1r3), Y2(:,c1r3), U1(:,c1r3),U2(:,c1r3)],'Ts',1,'TimeUnit','hours');

d4 = iddata([Y1(:,c2r1), Y2(:,c2r1)],[Y1(:,c2r1), Y2(:,c2r1), U1(:,c2r1),U2(:,c2r1)],'Ts',1,'TimeUnit','hours');
d5 = iddata([Y1(:,c2r2), Y2(:,c2r2)],[Y1(:,c2r2), Y2(:,c2r2), U1(:,c2r2),U2(:,c2r2)],'Ts',1,'TimeUnit','hours');
d6 = iddata([Y1(:,c2r3), Y2(:,c2r3)],[Y1(:,c2r3), Y2(:,c2r3), U1(:,c2r3),U2(:,c2r3)],'Ts',1,'TimeUnit','hours');

d7 = iddata([Y1(:,c3r1), Y2(:,c3r1)],[Y1(:,c3r1), Y2(:,c3r1), U1(:,c3r1),U2(:,c3r1)],'Ts',1,'TimeUnit','hours');
d8 = iddata([Y1(:,c3r2), Y2(:,c3r2)],[Y1(:,c3r2), Y2(:,c3r2), U1(:,c3r2),U2(:,c3r2)],'Ts',1,'TimeUnit','hours');
d9 = iddata([Y1(:,c3r3), Y2(:,c3r3)],[Y1(:,c3r3), Y2(:,c3r3), U1(:,c3r3),U2(:,c3r3)],'Ts',1,'TimeUnit','hours');

d10 = iddata([Y1(:,c4r1), Y2(:,c4r1)],[Y1(:,c4r1), Y2(:,c4r1), U1(:,c4r1),U2(:,c4r1)],'Ts',1,'TimeUnit','hours');
d11 = iddata([Y1(:,c4r2), Y2(:,c4r2)],[Y1(:,c4r2), Y2(:,c4r2), U1(:,c4r2),U2(:,c4r2)],'Ts',1,'TimeUnit','hours');
d12 = iddata([Y1(:,c4r3), Y2(:,c4r3)],[Y1(:,c4r3), Y2(:,c4r3), U1(:,c4r3),U2(:,c4r3)],'Ts',1,'TimeUnit','hours');

%totdatobj = merge(d1,d2,d3,d7,d8,d9); % 1 uM  arabinose/ 1 nM HSL netrecon set 
totdatobj = merge(d4,d5,d6,d10,d11,d12); % 1 mM arabinose / 1 uM HSL net recon set 
tm = t;

Ydim = 2;
