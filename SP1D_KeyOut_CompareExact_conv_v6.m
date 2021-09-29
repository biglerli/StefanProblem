function [l1vec,l2vec] =   SP1D_KeyOut_CompareExact_conv_v6(Mx,dt_init,T)
% Stefan Problem 1D
% exact solution is from RBCi paper with realistic coefficients 
% adding time adaptive 

% clear all
% clf(3)
% clf(4)
ploton = false;
printon = false;
%% Save Data
savedata = false;
if savedata 
    ploton = true;
    foldername = ['SP_GradSem21',datestr(now, 'ddmmmyyyyHHMM')];
    mkdir (foldername)
end
%% Save Movie
savemovie = false;
if savemovie 
    ploton = true;
    if ~savedata
        foldername = ['SP_GradSem21',datestr(now, 'ddmmmyyyyHHMM')];
        mkdir (foldername)
    end
    mov1 = VideoWriter([foldername '/SP1D.mov'],'MPEG-4');
    open(mov1)
end
ax=0;
bx=20;
% Mx=2000;

hx=(bx-ax)/Mx;
x=hx*(.5:Mx)+ax; x=x';


% T=200;
% dt_init=100*hx*hx;
dtmax=dt_init;


% uforv=cos(pi.*x);

Lw=306;
Lo=0;

ft_w = 0;
ft_o = 0;

scale =1; 


maxiter = 8;
acc = 1e-12;
tol = 1000;
largeiter = 0;
iter = 0;

% zero for material o and 1 for material w
 domvec = ones(size(x)); % domvec(find(x>0.2 & x<0.5))=0;
% domvec = 1-domvec;
Dom = domvec;
Lmat = Lw.*domvec+Lo.*(1-domvec);
ft_mat = ft_w.*domvec+ft_o.*(1-domvec);

u_av = [];
v_av = [];
phi_av = [];
% persistent bc
%%
example = 7;

bc = BoundCond(example,x,0);
uforv = bc.uforv;

% example = 1;
% 
% if example == 1
%     bc.left=1; bc.Ulbound=5;
%     bc.right=0;
% %     bc.Urbound=-1;
%     uforv=-5.*ones(size(x));
% %     if bc.left == 1
% %         uforv(1)=bc.Ulbound;
% %     end
% %     if bc.right == 1
% %         uforv(end)=bc.Urbound;
% %     end
% end
% 
% if example == 2
%     bc.left = 0;
%     bc.right = 0;
% end
% 
% if example == 3
%     bc.left=1;
%     bc.Ulbound=-1;
%     bc.right=1;
%     bc.Urbound=1;
%     uforv=-ones(size(x));
% end

c = cfun(uforv,domvec,ft_w,ft_o)./scale;
% v0=zeros(size(uforv));
%     for i =1:length(uforv)
%         v0(i) = c(i)*uforv(i)-max(0,min(uforv(i)-ft_mat(i),Lmat(i)/c(i)));%equation for 
%     end
% % v0=uforv;%cfun(uforv).*uforv+max(0,min(uforv,L));
v0=bc.H0;
v=v0;
u=uforv; % initial guess for u
un=u;
unk=u;
vn=v;
vnk=v;

dt=dt_init;
t_vec = 0;
t_count = 0;
t_cur=0;

bc = BoundCond(7,x,t_cur);
c = cfun(u,domvec,ft_w,ft_o)./scale;
k = kfun(u,domvec,ft_w,ft_o);
alpha = k./scale;
kx = kedgex(alpha, hx,bc,dt);
A = Amat(kx, Mx);
F=rhs(x,bc,Mx,dt,hx,kx,example,scale);
resV = zeros(size(vnk));
pdvnk = zeros(Mx,Mx);
% resplus = ones(2*Mx,maxiter);
% Allplus = 100.*ones(2*Mx,tsteps*maxiter);

%% initial plot
n=0;
% figure
%     plot(x,vn,'--','linewidth',3,'Color',[0.4660 0.6740 0.1880]);
% %       axis([a,b,-0.5,2]);
%     tt=dt*n;
%     title(['h = ',num2str(hx),', \tau = ',num2str(dt),', t = ',num2str(tt)],'fontsize',20);
% %     legend('enthalpy','fontsize',20);
%     xlabel('x','fontsize',12)
%     ylim([-3 6.5])
%     hold on
%     pause;
    
%%
k_up_vec = [];
c_up_vec = [];
u_avg_vec = [];
l1vec = [];
l2vec = [];
%% time looping
while t_cur<T

%     if iter == maxiter
%         disp(['max iter met, final time is ',num2str(n*dt)])
%         break
%     else
%         iter = 1;
%     end
    if (T-t_cur)<dt
        if T-t_cur<100*eps
            disp('done')
            break
        else
            dt = T-t_cur;
            t_cur = t_cur+dt;
%             kx = kedgex(alpha, hx, bc,dt);
%             A = Amat(kx, Mx);
%             F=rhs(x,bc,Mx,dt,hx,kx,example,scale);
        end
    elseif t_count>5
        disp(['max time count met, final time is ',num2str(t_cur)])
        break    
    elseif iter==maxiter 
        unk = un;
        vnk = vn;
        if t_count>1
            t_cur = t_cur-dt
        end
        dt=dt*.1;
        t_cur = t_cur+dt;
        t_count = t_count+1
        disp(['t_count=',num2str(t_count)])
        
%         c = cfun(un,domvec,ft_w,ft_o)./scale;
%         k = kfun(un,domvec,ft_w,ft_o);
%         alpha = k./scale;
%         kx = kedgex(alpha, hx, bc,dt);
%         A = Amat(kx, Mx);
%         F=rhs(x,bc,Mx,dt,hx,kx,example,scale);
    elseif iter<3         
        if dt<dtmax
            dt = dt*10
        end
        t_cur = t_cur+dt;
    end
    %%  set up matrix with new t_cur
    if printon
        disp(['tcur  = ',num2str(t_cur)])
    end
        bc = BoundCond(7,x,t_cur);
        c = cfun(un,domvec,ft_w,ft_o)./scale;
        k = kfun(un,domvec,ft_w,ft_o);

        alpha = k./scale;

        kx = kedgex(alpha, hx, bc,dt);
        A = Amat(kx, Mx);
        F = rhs(x,bc,Mx,dt,hx,kx,example,scale);
    %%    
    iter = 1;
    if printon
        disp(['time step is ',num2str(t_cur)])
    end
    while iter<maxiter
        if printon
            disp(['iter=',num2str(iter)])
        end
        % TODO build A matrix based on K(u)
%         c = cfun(un,domvec,ft_w,ft_o);
%         k = kfun(un,domvec,ft_w,ft_o);
%         
%         kx = kedgex(k, hx, bc,dt);
%         A = Amat(kx, Mx);
%         F=rhs(x,bc,Mx,dt,hx,kx,example,scale);
        
        %full vector
        Allnk = [unk; vnk];
        Allplus(:,iter) = Allnk;
        
        %residual
        resSP = hx*(vnk-vn)+A*unk-F; % SP multiplied by hx*dt
        for i =1:length(vnk)
            resV(i) = vnk(i)-c(i)*unk(i)-max(0,min(vnk(i)-c(i)*ft_mat(i),Lmat(i)/scale));%equation for 
        end
        res = [resSP; resV];
%         resplus(:,iter) = res;
        %jacobian
        
        
        
        if norm(res,inf)<tol
            pdunk = -diag(c);
%             pdunk = -eye(Mx);
            for i=1:length(vnk)
                if vnk(i)-c(i)*ft_mat(i)<=Lmat(i)/scale && vnk(i)-c(i)*ft_mat(i)>=0
                    pdvnk(i,i) = 0;
                else
                    pdvnk(i,i) = 1;
                end
            end
                
            pdvnk = sparse(pdvnk);
                        
            %jacobi: SP/U, SP/V; V/V V/U
            Jac =[[A;pdunk],[sparse(hx.*eye(Mx));pdvnk]];
            
            s= Jac\res;
            Allnk = Allnk-s;
            unk  = Allnk(1:Mx);
            vnk  = Allnk(Mx+1:end);
            
            
           
            if norm(res,inf)<acc
                largeiter = max(largeiter,iter);
%                 t_cur = t_cur+dt;
                t_vec = [t_vec;dt];
                if iter<3 && (T-t_cur)>10*eps
                    if dt<dtmax
                        dt = dt*10
                    end
                end
                t_count = 0;
                All = Allnk;
                un = All(1:Mx);
                vn = All(Mx+1:end);
                
                [l1t,l2t] = exactSolErrL2(ax,bx,x,hx,un,t_cur);%exactSolErr
                l1vec = [l1vec;l1t];
                l2vec = [l2vec;l2t];
                
                %% for upscale
%                 u_av = [u_av; mean(un)];
%                 v_av = [v_av; mean(vn)];
%                 phi_av = [phi_av;mean(sign(un))];
                
                %% for non-time dependent
%                 c = cfun(un,domvec,ft_w,ft_o)./scale;
%                 k = kfun(un,domvec,ft_w,ft_o);
%                 
%                 alpha = k./scale;
% 
%                 kx = kedgex(alpha, hx, bc,dt);
%                 A = Amat(kx, Mx);
%                 F = rhs(x,bc,Mx,dt,hx,kx,example,scale);
                
                %% for upscale
%                 if mod(length(t_vec)-1,1) == 0 
%                     k_up_vec=[k_up_vec;upscale_k(k,x,Mx,hx,bx-ax,t_cur)];
%                     c_up_vec=[c_up_vec;upscale_c(c,x,Mx,hx,bx-ax,t_cur)];
%                     u_avg_vec = [u_avg_vec;mean(un)];
%                 end
                plotint = 67000;
                if mod(floor(t_cur+eps),plotint) == 0 && ploton
%                 if mod(length(t_vec)-2,6000) == 0 && ploton 
%                 if abs(t_cur-0.48)<100*eps
                    figure(2)
%                     clf(3)
                    hold on
%                     plot(x,vn,'-.','linewidth',3,'Color',[0 0.4470 0.7410]);
                    plot(x,un,'o','linewidth',3); %,'Color',[0.8500 0.3250 0.0980]);
                %       axis([a,b,-0.5,2]);
                    tt=t_cur/(3600);
%                     format bank
%                     title([' t = ',num2str(tt,'%.1f'),'[hours]'],'fontsize',20);
%                     legend('enthalpy','temperature','fontsize',20);
%                     legend('v_0','v','u','fontsize',20);
                    xlabel('x','fontsize',20)
                    set(gca,'FontSize',20)
                    
                    figure(1)
%                     clf(4)
                    hold on
                    plot(x,vn,'x','linewidth',3); %,'Color',[0 0.4470 0.7410]);
%                     plot(x,un,'-','linewidth',3,'Color',[0.8500 0.3250 0.0980]);
                %       axis([a,b,-0.5,2]);
                    tt=t_cur/(3600);
%                     format bank
%                     title([' t = ',num2str(tt,'%.1f'),'[hours]'],'fontsize',20);
%                     legend('enthalpy','temperature','fontsize',20);
%                     legend('v_0','v','u','fontsize',20);
                    xlabel('x','fontsize',20)
                    set(gca,'FontSize',20)

%                     ylim([-10 400])
                    hold off
%                     pause;
                   if savemovie 
                       frame = getframe(3);
                       writeVideo(mov1,frame);
                   end
                end
                break
            else
                iter = iter+1;
%                 disp(k(6))
%                 disp(unk(1))
            end
                
        else
            error('initial guess too far from zero')
        end
    end
end
    
q_1 = kx(1)*(un(1)-bc.Ulbound)*(bx-ax);

if iter == maxiter
    disp('maxiter :/')
end
if ploton
    figure
    plot(x,vn,'--','linewidth',3);
    figure
    plot(x,un,':','linewidth',3);
    %       axis([a,b,-0.5,2]);
    tt=t_cur;
    % title(sprintf('t=%g',tt),'fontsize',20);
    % legend('enthalpy','temperature','fontsize',20);
    legend('v','u','v_{final}','u_{final}','fontsize',20);
    % xlabel('x','fontsize',12)
    ylim([-3 6.5])

    figure
    plot(u_av,v_av,'linewidth',3)
    xlabel('u','fontsize',15)
    ylabel('v','fontsize',15)
end

%% error
        c_i = 1.90;
        c_w = 4.19;
        L = 306;
        k_i = 0.023;
        k_w = 0.0058;

        d_i = k_i/c_i;
        d_w = k_w/c_w;

        v_s = -5e-5;

        B = -594;
        alpha = v_s/d_w;
        beta = v_s/d_i;
      s = 15+v_s*t_cur;
            Hend = x.*0;
            uend = x.*0;
            for i=1:length(x)
                % above freezing
                if x(i)<=s 
                    Hend(i)=-B+(B+L)*exp(alpha*(v_s*t_cur-x(i)+15));
                    uend(i)=(Hend(i)-L)/c_w;
                % below freezing
                else
                    Hend(i)=-B+(B)*exp(beta*(v_s*t_cur-x(i)+15));
                    uend(i)=(Hend(i))/c_i;
                end
            end
%        l1norm = norm(un-uend,1)*hx;
%        l2norm = norm(un-uend,2)*sqrt(hx)
%% movie
if savemovie
    close(mov1)
end
%% save data
if savedata
    datastr = ["hx";"hy";"T";"largeiter"];
    datanum = [hx;hy;T;largeiter];
    fileID = fopen([foldername '/RunDetail.txt'],'w');
    for i=1:length(datastr)
        fprintf(fileID,'%s %f\n',datastr(i),datanum(i));
    end
    fclose(fileID);
end

save('../SPuavgvec5.mat','u_avg_vec')
save('../SPkupvec5.mat','k_up_vec')
save('../SPcupvec5.mat','c_up_vec')
% function c=cfun(u)
%     c=zeros(size(u));
%     for i=1:length(u)
%         if u(i) <=0
%             c(i)=1;
%         elseif u(i)<100
%             c(i)=1;
%         else
%             error('water is boiling')
%         end
%     end
% end
% 
% function k=kfun(u)
%     k=zeros(size(u));
%     for i = 1:length(u)
%         if u(i) <= 0
%             k(i) = 2;
%         elseif u(i)<100
%             k(i) = 1;
%         else
%             error('water is boiling')
%         end
%     end
% end
end 
function [l1,l2] = exactSolErrL2(ax,bx,x,hx,un,t_cur)
%% error
        hx=hx/10;
        Mx=length(x).*10;
%         disp(['mx=',num2str(Mx)])
        
        x_e=hx*(.5:Mx)+ax; x_e=x_e';

        c_i = 1.90;
        c_w = 4.19;
        L = 306;
        k_i = 0.023;
        k_w = 0.0058;

        d_i = k_i/c_i;
        d_w = k_w/c_w;

        v_s = -5e-5;

        B = -594;
        alpha = v_s/d_w;
        beta = v_s/d_i;
      s = 15+v_s*t_cur;
            Hend = x_e.*0;
            uend = x_e.*0;
            for i=1:length(x_e)
                % above freezing
                if x_e(i)<=s 
                    Hend(i)=-B+(B+L)*exp(alpha*(v_s*t_cur-x_e(i)+15));
                    uend(i)=(Hend(i)-L)/c_w;
                % below freezing
                else
                    Hend(i)=-B+(B)*exp(beta*(v_s*t_cur-x_e(i)+15));
                    uend(i)=(Hend(i))/c_i;
                end
            end
            
            err_vec = [];
            
%             for j=1:length(uend)
%                 disp(['j=',num2str(j)])
                for m=1:length(x)
                    for k=1:10
                        err_vec(10*(m-1)+k) = uend(10*(m-1)+k)-un(m);
%                         disp(uend(10*(m-1)+k))
%                         j=j+1;
                    end
%                     j=j-1;
                end 
%             end
       l1 = norm(err_vec,1)*hx;
       l2 = norm(err_vec,2)*sqrt(hx);
end
function [l1,l2] = exactSolErr(ax,bx,x,hx,un,t_cur)
%% error
        c_i = 1.90;
        c_w = 4.19;
        L = 306;
        k_i = 0.023;
        k_w = 0.0058;

        d_i = k_i/c_i;
        d_w = k_w/c_w;

        v_s = -5e-5;

        B = -594;
        alpha = v_s/d_w;
        beta = v_s/d_i;
      s = 15+v_s*t_cur;
            Hend = x.*0;
            uend = x.*0;
            for i=1:length(x)
                % above freezing
                if x(i)<=s 
                    Hend(i)=-B+(B+L)*exp(alpha*(v_s*t_cur-x(i)+15));
                    uend(i)=(Hend(i)-L)/c_w;
                % below freezing
                else
                    Hend(i)=-B+(B)*exp(beta*(v_s*t_cur-x(i)+15));
                    uend(i)=(Hend(i))/c_i;
                end
            end
       l1 = norm(un-uend,1)*hx;
       l2 = norm(un-uend,2)*sqrt(hx);
end
function bc = BoundCond(example,x,t)
    if example == 1
        bc.left=1; bc.Ulbound=20;
        bc.right=0;
    %     bc.Urbound=-1;
        bc.uforv=-10.*ones(size(x));
    %     if bc.left == 1
    %         uforv(1)=bc.Ulbound;
    %     end
    %     if bc.right == 1
    %         uforv(end)=bc.Urbound;
    %     end
    end

    if example == 2
        bc.left = 0;
        bc.right = 0;
    end

    if example == 3
        bc.left=1;
        bc.Ulbound=-1;
        bc.right=1;
        bc.Urbound=1;
        bc.uforv=-ones(size(x));
    end
    
    if example == 4
        % for upscaling
        bc.left=1; bc.Ulbound=1;
        bc.right=1; bc.Urbound=0;
    %     bc.Urbound=-1;
        bc.uforv=-ones(size(x));
    end
    
    if example == 5
        bc.left=1; bc.Ulbound=5;
        bc.right=1;bc.Ulbound=5;
    %     bc.Urbound=-1;
        bc.uforv=-ones(size(x));
    %     if bc.left == 1
    %         uforv(1)=bc.Ulbound;
    %     end
    %     if bc.right == 1
    %         uforv(end)=bc.Urbound;
    %     end
    end
    
        %for comparison
    if example == 6
        % for upscaling
        bc.left=1; bc.Ulbound=10;
        bc.right=0; 
    %     bc.Urbound=-1;
        bc.uforv=-5.*ones(size(x));
    end
    
        c_i = 1.90;
        c_w = 4.19;
        L = 306;
        k_i = 0.023;
        k_w = 0.0058;

        d_i = k_i/c_i;
        d_w = k_w/c_w;

        v_s = -5e-5;

        B = -594;
        alpha = v_s/d_w;
        beta = v_s/d_i;
        
    if example == 7
        % for upscaling
        bc.left=1;  
        bc.right=1; 
                Hl=-B+(B+L)*exp(alpha*(v_s*t+15));
                bc.Ulbound=(Hl-L)/c_w;
            % below freezing
                Hr=-B+(B)*exp(beta*(v_s*t-5));
                bc.Urbound=(Hr)/c_i;
        %     bc.Urbound=-1;
        if t==0
            H0 = x.*0;
            u0 = x.*0;
            s=15;
            for i=1:length(x)
                % above freezing
                if x(i)<=s 
                    H0(i)=-B+(B+L)*exp(alpha*(v_s*t-x(i)+15));
                    u0(i)=(H0(i)-L)/c_w;
                % below freezing
                else
                    H0(i)=-B+(B)*exp(beta*(v_s*t-x(i)+15));
                    u0(i)=(H0(i))/c_i;
                end
            end
            bc.uforv=u0;
            bc.H0=H0;
        end
    end
end

function k_up = upscale_k(k,x,Mx,hx,x_side,t)
    bcup = BoundCond(4,x,t);
    bcup.Ulbound = bcup.Ulbound.*x_side;
%     c = cfun(u,domvec,ft_w,ft_o)./scale;
%     k = kfun(u,domvec,ft_w,ft_o);
    kx = kedgex(k,hx,bcup,1);
    A_up = Amat(kx, Mx);
    F_up=rhs(x,bcup,Mx,1,hx,kx,4,1);
    
%     k = kfun(u,myin,Dom,ft_w,ft_o);
%     kx = kedgex(k,hx,hy,Mx,My,1,BC);
%     ky = kedgey(k,hx,hy,Mx,My,1,BC);
%     A_x = Amat(kx,ky,BC,hx,hy,Mx,My);
%     Fc = zeros(Mx*My,1); %forcing function
%     F = rhs(x,y,BC,Mx,My,hx,hy,kx,ky,1,Fc);

    up = A_up\F_up;

%     % unpack the solution
%     for i=1:Mx
%         for j=1:My
%             me=myin(i,j);
%             U_temp_p(i,j)=u_temp(me);
%         end 
%     end

%     k_av_harm_x = harmmean(k)
%     k_av_x = mean(k)
    % q_1 = kx(1,:)*(U_temp_p(1,:)-BC.Ulbound);
    % q_end = kx(end,:)*(BC.Urbound-U_temp_p(end,:));

    q_in_x = -kx(1)*(up(1)-bcup.Ulbound);
%     q_out_x = -kx(end)*(bcup.Urbound-up(end));
    
    k_up = q_in_x;

    
end

function c_up = upscale_c(c,x,Mx,hx,x_side,t)
    bcup = BoundCond(4,x,t);
    bcup.Ulbound = bcup.Ulbound.*x_side;

    c_up = sum(c)*hx/x_side;

    
end

function c=cfun(u,domvec,ft_w,ft_o)
    c=zeros(size(u));
%     cw_s =2.1; %ice
%     cw_l =2.2; %water
% %     cw_l =4.2; %water
%     co_s =3; %other solid
%     co_l =3.2; %other liquid
    
    cw_s =1.9; %ice
    cw_l =4.19; %water
    co_l =0.7; %water
    co_s =0.7; %other solid
%     co_l =3.2; %other liquid
    
    for i=1:length(domvec)
        if domvec(i)==1
            if u(i) <=ft_w
                c(i)=cw_s;
            elseif u(i)<100
                c(i)=cw_l;
            else
                error('water is boiling')
            end
        else
            if u(i) <=ft_o
                c(i)=co_s;
            elseif u(i)<100
                c(i)=co_l;
            else
                error('other is boiling')
            end
        end
    end
    
%     c = c./scale;
end

function k = kfun(u,domvec,ft_w,ft_o)
     k=zeros(size(u));
%     kw_s =2.2; %ice
% %     kw_l =2.2; %water
%     kw_l =0.6; %water
%     kw_s =1.6; %ice
% %     kw_l =2.2; %water
%     kw_l =1.8; %water
%     ko_s =3; %other solid
%     ko_l =1; %other liquid
    
    
    kw_s =.023; %ice
    kw_l = .0058; %water
    ko_s =1.4; %other solid
    ko_l =1.4; %other liquid
    
    
    for i=1:length(domvec)
        if domvec(i)==1
            if u(i) <=ft_w
                k(i)=kw_s;
            elseif u(i)<100
                k(i)=kw_l;
            else
                error('other is boiling')
            end
%             if u(i) <=ft_w-eps
%                 k(i)=kw_s;
%             elseif u(i) >=ft_w+eps
%                 k(i)=kw_l;
%             else
%                 k(i)=2/(1/kw_l+1/kw_s);
%             end
        else
            if u(i) <=ft_o
                k(i)=ko_s;
            elseif u(i)<100
                k(i)=ko_l;
            else
                error('other is boiling')
            end
        end
    end

end

function F=Ffun(x,example)
    if example == 1
        F = zeros(size(x));
%         F = ones(size(x));
    end
    if example == 2
        F = zeros(size(x));
%         F = 10*(ones(size(x))-2*around(x));
    end
        
    if example == 3
        F = zeros(size(x));
        F(find(x>0.7 & x<0.8))=10;
    end 
    
    if example == 4
        F = zeros(size(x));
    end
    
    if example == 6
        F = zeros(size(x));
    end
    
    if example == 7
        F = 0.*ones(size(x));
    end
end

function kx=kedgex(k,hx,bc,dt)
    kx = zeros(length(k)+1,1);
    for i =2:length(k)
        if k(i-1)*k(i)==0
            kx(i)=0;
        else
            kx(i)=2*dt/(hx/k(i-1)+hx/k(i));
        end
        
        if bc.left == 1
            kx(1) = 2*dt*k(1)/hx;
        else
            kx(1)=0;
        end
        
        if bc.right == 1
            kx(end) = 2*dt*k(end)/hx;
        else
            kx(end)=0;
        end
    end
      
end

function A = Amat(kx,Mx)
    A=zeros(Mx,Mx);
    for l=2:Mx
        cL = l-1;
        cR = l;
        A(cL,cL) = A(cL,cL)+kx(l);
        A(cR,cR) = A(cR,cR)+kx(l);
        A(cR,cL) = A(cR,cL)-kx(l);
        A(cL,cR) = A(cL,cR)-kx(l);
    end
    A(1,1)=A(1,1)+kx(1);
    A(end,end)=A(end,end)+kx(end);
    A = sparse(A);
end

function frhs = rhs(x,bc,Mx,dt,hx,kx,example,scale)
    fc = Ffun(x,example);
    frhs = (hx*dt/scale).*fc;
    if bc.left ==1
        frhs(1)=frhs(1)+bc.Ulbound*kx(1);
    end
    if bc.right==1
        frhs(end)=frhs(end)+bc.Urbound*kx(end);
    end   
end

