%%
tstep=0.005;
nstep=10/tstep;
%%
mi=[1,9,6,4,1,0.6,0.5];
Jip(:,:,1)=diag([1,1,1]);
Jip(:,:,2)=diag([0.01,0.01,0.02]);
Jip(:,:,3)=diag([0.05,0.06,0.01]);
Jip(:,:,4)=diag([0.4,0.4,0.01]);
Jip(:,:,5)=diag([0.001,0.0005,0.001]);
Jip(:,:,6)=diag([0.0005,0.0005,0.0002]);
Jip(:,:,7)=diag([0.003,0.001,0.002]);
%
qij(:,1,2)=[0;0];
qij(:,2,3)=[0;0];
qij(:,3,4)=[0;0];
qij(:,4,5)=[0;0];
qij(:,5,6)=[0;0];
qij(:,6,7)=[0;0];
%dqij=[dtheta;dalpha]
dqij(:,1,2)=[0;0];
dqij(:,2,3)=[0;0];
dqij(:,3,4)=[0;0];
dqij(:,4,5)=[0;0];
dqij(:,5,6)=[0;0];
dqij(:,6,7)=[0;0];
%
ddqij(:,1,2)=[0;0];
ddqij(:,2,3)=[0;0];
ddqij(:,3,4)=[0;0];
ddqij(:,4,5)=[0;0];
ddqij(:,5,6)=[0;0];
ddqij(:,6,7)=[0;0];
%%
%sijp & sjip
sijp(:,1,2)=[0;0;0];sjip(:,2,1)=[0;0;0];
sijp(:,2,3)=[0.1;0;0.1];sjip(:,3,2)=[0;0;0];
sijp(:,3,4)=[0;0;0];sjip(:,4,3)=[0;0;0];
sijp(:,4,5)=[0;0.3;0];sjip(:,5,4)=[0;-0.3;0];
sijp(:,5,6)=[0;0.1;0];sjip(:,6,5)=[0;-0.06;0];
sijp(:,6,7)=[0;0.04;0];sjip(:,7,6)=[0;-0.1;0];
%Cij & Cji
Cij(:,:,1,2)=eye(3);Cji(:,:,2,1)=eye(3);
Cij(:,:,2,3)=eye(3);Cji(:,:,3,2)=eye(3);
Cij(:,:,3,4)=eye(3);Cji(:,:,4,3)=eye(3);
Cij(:,:,4,5)=eye(3);Cji(:,:,5,4)=eye(3);
Cij(:,:,5,6)=eye(3);Cji(:,:,6,5)=eye(3);
Cij(:,:,6,7)=eye(3);Cji(:,:,7,6)=eye(3);
%dijpp & djipp
dijpp(:,1,2)=[0;0;1];djipp(:,2,1)=[0;0;1];
dijpp(:,2,3)=[1;0;0];djipp(:,3,2)=[1;0;0];
dijpp(:,3,4)=[0;1;0];djipp(:,4,3)=[0;1;0];
dijpp(:,4,5)=[0;1;0];djipp(:,5,4)=[0;1;0];
dijpp(:,5,6)=[0;1;0];djipp(:,6,5)=[0;1;0];
dijpp(:,6,7)=[0;1;0];djipp(:,7,6)=[0;1;0];
%
hijpp(:,1,2)=[0;0;1];hjipp(:,2,1)=[0;0;1];
hijpp(:,2,3)=[1;0;0];hjipp(:,3,2)=[1;0;0];
hijpp(:,3,4)=[0;1;0];hjipp(:,4,3)=[0;1;0];
hijpp(:,4,5)=[0;1;0];hjipp(:,5,4)=[0;1;0];
hijpp(:,5,6)=[0;0;1];hjipp(:,6,5)=[0;0;1];
hijpp(:,6,7)=[0;1;0];hjipp(:,7,6)=[0;1;0];
%%
%Aijpp
Aijpp(:,:,1,2)=CgaCal(0,0,qij(1,1,2));
Aijpp(:,:,2,3)=CgaCal(qij(1,2,3),0,0);
Aijpp(:,:,3,4)=CgaCal(0,0,0);
Aijpp(:,:,4,5)=CgaCal(0,qij(1,4,5),0);
Aijpp(:,:,5,6)=CgaCal(0,0,qij(1,5,6));
Aijpp(:,:,6,7)=CgaCal(0,qij(1,6,7),0);
%Ai
Ai(:,:,1)=eye(3);
for i = 1:6
    Ai(:,:,i+1)=Ai(:,:,i)*Cij(:,:,i,i+1)*Aijpp(:,:,i,i+1)*Cji(:,:,i+1,i);
end
% Ai(:,:,2)=Ai(:,:,1)*Cij(:,:,1,2)*Aijpp(:,:,1,2)*Cji(:,:,2,1);
% Ai(:,:,3)=Ai(:,:,2)*Cij(:,:,2,3)*Aijpp(:,:,2,3)*Cji(:,:,3,2);
% Ai(:,:,4)=Ai(:,:,3)*Cij(:,:,3,4)*Aijpp(:,:,3,4)*Cji(:,:,4,3);
% Ai(:,:,5)=Ai(:,:,4)*Cij(:,:,4,5)*Aijpp(:,:,4,5)*Cji(:,:,5,4);
% Ai(:,:,6)=Ai(:,:,5)*Cij(:,:,5,6)*Aijpp(:,:,5,6)*Cji(:,:,6,5);
% Ai(:,:,7)=Ai(:,:,6)*Cij(:,:,6,7)*Aijpp(:,:,6,7)*Cji(:,:,7,6);
%
for i = 1:7
    Ji(:,:,i)=Ai(:,:,i)*Jip(:,:,i)*Ai(:,:,i)';
    Mi(:,:,i)=[diag([mi(i),mi(i),mi(i)]),zeros(3);zeros(3),Ji(:,:,i)];
end
%
for i = 1:6
    dij(:,i,i+1)=Ai(:,:,i)*Cij(:,:,i,i+1)*dijpp(:,i,i+1);
    hij(:,i,i+1)=Ai(:,:,i)*Cij(:,:,i,i+1)*hijpp(:,i,i+1);
    sij(:,i,i+1)=Ai(:,:,i)*sijp(:,i,i+1);
    sji(:,i+1,i)=Ai(:,:,i)*Cij(:,:,i,i+1)*Aijpp(:,:,i,i+1)*Cji(:,:,i+1,i)'*sjip(:,i+1,i);
end
%
wi(:,1)=[0;0;0];wimat(:,:,1)=Vec2Mat(wi(:,1));
for i = 2:7
    wi(:,i)=wi(:,i-1)+hij(:,i-1,i)*dqij(1,i-1,i);
    wimat(:,:,i)=Vec2Mat(wi(:,i));
end
for i = 1:6
    dsij(:,i,i+1)=wimat(:,:,i)*sij(:,i,i+1);
    dsji(:,i+1,i)=wimat(:,:,i+1)*sji(:,i+1,i);
    ddij(:,i,i+1)=wimat(:,:,i)*dij(:,i,i+1);
    dhij(:,i,i+1)=wimat(:,:,i)*hij(:,i,i+1);
end
%rij
for i = 1:6
    rij(:,i,i+1)=sij(:,i,i+1)+dij(:,i,i+1)*qij(2,i,i+1)-sji(:,i+1,i);
    drij(:,i,i+1)=dsij(:,i,i+1)+ddij(:,i,i+1)*qij(2,i,i+1)+dij(:,i,i+1)*dqij(2,i,i+1)-dsji(:,i+1,i);
end
%ri
ri(:,1)=[0;0;0];
for i = 2:7
    ri(:,i)=ri(:,i-1)+rij(:,i-1,i);
end
%
for i = 1:6
    rijmat(:,:,i,i+1)=Vec2Mat(rij(:,i,i+1));
    hijmat(:,:,i,i+1)=Vec2Mat(hij(:,i,i+1));

    drijmat(:,:,i,i+1)=Vec2Mat(drij(:,i,i+1));
    dhijmat(:,:,i,i+1)=Vec2Mat(dhij(:,i,i+1));
end
% for i = 1:6
%     Bij1(:,:,i,i+1)=[eye(3),-rijmat(:,:,i,i+1);zeros(3),eye(3)];
%     Bij2(:,:,i,i+1)=[-hijmat(:,:,i,i+1)*sji(:,i+1,i),dij(:,i,i+1);hij(:,i,i+1),zeros(3,1)];
% 
%     dBij1(:,:,i,i+1)=[zeros(3),-drijmat(:,:,i,i+1);zeros(3),zeros(3)];
%     dBij2(:,:,i,i+1)=[-dhijmat(:,:,i,i+1)*sji(:,i+1,i)-hijmat(:,:,i,i+1)*dsji(:,i+1,i),ddij(:,i,i+1);dhij(:,i,i+1),zeros(3,1)];
% end
%
Yi(:,1)=[0,0,0,0,0,0]';dYi(:,1)=[0,0,0,0,0,0]';
% for i = 1:6
%     Yi(:,i+1)=Bij1(:,:,i,i+1)*Yi(:,i)+Bij2(:,:,i,i+1)*dqij(:,i,i+1);
%     Dij(:,i,i+1)=dBij1(:,:,i,i+1)*Yi(:,i)+dBij2(:,:,i,i+1)*dqij(:,i,i+1);
%     dYi(:,i+1)=Bij1(:,:,i,i+1)*dYi(:,i)+Bij2(:,:,i,i+1)*ddqij(:,i,i+1)+Dij(:,i,i+1);
% end
%%
%%
euleranglei(:,1)=[0;0;0];
t=0;
for step=1:nstep
    t=step*tstep;
    fi(:,1)=[0;0;0;0;0;0];
    fi(:,2)=[0;0;0;0;0;0.12*sin(pi*t/5)];
    fi(:,3)=[0;0;0;0;3*t/5+14.5;0];
    fi(:,4)=[0;0;15*cos(pi*t/10)-15;0;0;0];
    fi(:,5)=[0;0;0;0;0;0.4*cos(pi/10*t)-0.4];
    fi(:,6)=[0;0;0;0;0.2*cos(pi*t/10)+1.1;0];
    fi(:,7)=[0;0;0;0;0;0.00022*sin(pi/5*(t-5))+0.00004];

%     fi(:,1)=[0;0;0;0;0;0];
%     fi(:,2)=[0;0;0;0;0;0.12*sin(pi*t/5)];
%     fi(:,3)=[0;0;0;3*t/5+14.5;0;0];
    fi(:,4)=[0;15*cos(pi*t/10)-15;0;0;0;0];
%     fi(:,5)=[0;0;0;0;0.4*cos(pi/10*t)-0.4;0];
%     fi(:,6)=[0;0;0;0;0;0.2*cos(pi*t/10)+1.1];
%     fi(:,7)=[0;0;0;0;0.00022*sin(pi/5*(t-5))-0.00004;0];

    fi=zeros(6,7);
    fi(:,4)=[0;15*cos(pi*t/5)-15;0;0;0;0];
    fi(:,4)=[0;6.1;0;0;0;0];
    
%Aijpp
    Aijpp(:,:,1,2)=CgaCal(0,0,qij(1,1,2));
    Aijpp(:,:,2,3)=CgaCal(qij(1,2,3),0,0);
    Aijpp(:,:,3,4)=CgaCal(0,0,0);
    Aijpp(:,:,4,5)=CgaCal(0,qij(1,4,5),0);
    Aijpp(:,:,5,6)=CgaCal(0,0,qij(1,5,6));
    Aijpp(:,:,6,7)=CgaCal(0,qij(1,6,7),0);
    for i = 1:6
        Ai(:,:,i+1)=Ai(:,:,i)*Cij(:,:,i,i+1)*Aijpp(:,:,i,i+1)*Cji(:,:,i+1,i);
        
        dij(:,i,i+1)=Ai(:,:,i)*Cij(:,:,i,i+1)*dijpp(:,i,i+1);
        hij(:,i,i+1)=Ai(:,:,i)*Cij(:,:,i,i+1)*hijpp(:,i,i+1);
        sij(:,i,i+1)=Ai(:,:,i)*sijp(:,i,i+1);
        sji(:,i+1,i)=Ai(:,:,i)*Cij(:,:,i,i+1)*Aijpp(:,:,i,i+1)*Cji(:,:,i+1,i)'*sjip(:,i+1,i);
        
        wi(:,i+1)=wi(:,i)+hij(:,i,i+1)*dqij(1,i,i+1);
        wimat(:,:,i+1)=Vec2Mat(wi(:,i+1));
        hijmat(:,:,i,i+1)=Vec2Mat(hij(:,i,i+1));
        
        dsij(:,i,i+1)=wimat(:,:,i)*sij(:,i,i+1);
        dsji(:,i+1,i)=wimat(:,:,i)*sji(:,i+1,i)+hijmat(:,:,i,i+1)*sji(:,i+1,i)*dqij(1,i,i+1);
        ddij(:,i,i+1)=wimat(:,:,i)*dij(:,i,i+1);
        dhij(:,i,i+1)=wimat(:,:,i)*hij(:,i,i+1);

        rij(:,i,i+1)=sij(:,i,i+1)+dij(:,i,i+1)*qij(2,i,i+1)-sji(:,i+1,i);
        drij(:,i,i+1)=dsij(:,i,i+1)+ddij(:,i,i+1)*qij(2,i,i+1)+dij(:,i,i+1)*dqij(2,i,i+1)-dsji(:,i+1,i);
        ri(:,i+1)=ri(:,i)+rij(:,i,i+1);

        rijmat(:,:,i,i+1)=Vec2Mat(rij(:,i,i+1));
        hijmat(:,:,i,i+1)=Vec2Mat(hij(:,i,i+1));

        drijmat(:,:,i,i+1)=Vec2Mat(drij(:,i,i+1));
        dhijmat(:,:,i,i+1)=Vec2Mat(dhij(:,i,i+1));

        Bij1(:,:,i,i+1)=[eye(3),-rijmat(:,:,i,i+1);zeros(3),eye(3)];
%         Bij2(:,:,i,i+1)=[-hijmat(:,:,i,i+1)*sji(:,i+1,i),dij(:,i,i+1);hij(:,i,i+1),zeros(3,1)];
% 
        dBij1(:,:,i,i+1)=[zeros(3),-drijmat(:,:,i,i+1);zeros(3),zeros(3)];
%         dBij2(:,:,i,i+1)=[-dhijmat(:,:,i,i+1)*sji(:,i+1,i)-hijmat(:,:,i,i+1)*dsji(:,i+1,i),ddij(:,i,i+1);dhij(:,i,i+1),zeros(3,1)];

        if i==3
            Bij2(:,i,i+1)=[dij(:,i,i+1);zeros(3,1)];
            dBij2(:,i,i+1)=[ddij(:,i,i+1);zeros(3,1)];
            Yi(:,i+1)=Bij1(:,:,i,i+1)*Yi(:,i)+Bij2(:,i,i+1)*dqij(2,i,i+1);
            Dij(:,i,i+1)=dBij1(:,:,i,i+1)*Yi(:,i)+dBij2(:,i,i+1)*dqij(2,i,i+1);
        else
            Bij2(:,i,i+1)=[-hijmat(:,:,i,i+1)*sji(:,i+1,i);hij(:,i,i+1)];
            dBij2(:,i,i+1)=[-dhijmat(:,:,i,i+1)*sji(:,i+1,i)-hijmat(:,:,i,i+1)*dsji(:,i+1,i);dhij(:,i,i+1)];
            Yi(:,i+1)=Bij1(:,:,i,i+1)*Yi(:,i)+Bij2(:,i,i+1)*dqij(1,i,i+1);
            Dij(:,i,i+1)=dBij1(:,:,i,i+1)*Yi(:,i)+dBij2(:,i,i+1)*dqij(1,i,i+1);
        end
%         Yi(:,i+1)=Bij1(:,:,i,i+1)*Yi(:,i)+Bij2(:,:,i,i+1)*dqij(:,i,i+1);
%         Dij(:,i,i+1)=dBij1(:,:,i,i+1)*Yi(:,i)+dBij2(:,:,i,i+1)*dqij(:,i,i+1);
        
        %dYi(:,i+1)=Bij1(:,:,i,i+1)*dYi(:,i)+Bij2(:,:,i,i+1)*ddqij(:,i,i+1);
    end
    
    for i = 1:7
        Ji(:,:,i)=Ai(:,:,i)*Jip(:,:,i)*Ai(:,:,i)';
        Mi(:,:,i)=[diag([mi(i),mi(i),mi(i)]),zeros(3);zeros(3),Ji(:,:,i)];
    end
    
    for i = 1:7
        Qi(:,i)=[Ai(:,:,i),zeros(3);zeros(3),Ai(:,:,i)]*fi(:,i)+[0;0;0;-wimat(:,:,i)*Ji(:,:,i)*wi(:,i)];
    end
    
    Ki=zeros(6,6,7);Li=zeros(6,7);
%     for i = 7:-1:2%反向迭代语法
%         Ki(:,:,i-1)=Bij1(:,:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij1(:,:,i-1,i)-...
%             Bij1(:,:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij2(:,:,i-1,i)*...
%             (Bij2(:,:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij2(:,:,i-1,i))^(-1)*...
%             (Bij2(:,:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij1(:,:,i-1,i));
% 
%         Li(:,i-1)=Bij1(:,:,i-1,i)'*(-(Mi(:,:,i)+Ki(:,:,i))*Dij(:,i-1,i)+(Qi(:,i)+Li(:,i))+...
%             (Mi(:,:,i)+Ki(:,:,i))*Bij2(:,:,i-1,i)*(Bij2(:,:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij2(:,:,i-1,i))^(-1)*...
%             (Bij2(:,:,i-1,i)'*((Mi(:,:,i)+Ki(:,:,i))*Dij(:,i-1,i)-(Qi(:,i)+Li(:,i)))));
%         
%     end

    for i = 7:-1:2%反向迭代语法
        Ki(:,:,i-1)=Bij1(:,:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij1(:,:,i-1,i)-...
            Bij1(:,:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij2(:,i-1,i)*...
            (Bij2(:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij2(:,i-1,i))^(-1)*...
            (Bij2(:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij1(:,:,i-1,i));

        Li(:,i-1)=Bij1(:,:,i-1,i)'*(-(Mi(:,:,i)+Ki(:,:,i))*Dij(:,i-1,i)+(Qi(:,i)+Li(:,i))+...
            (Mi(:,:,i)+Ki(:,:,i))*Bij2(:,i-1,i)*(Bij2(:,i-1,i)'*(Mi(:,:,i)+Ki(:,:,i))*Bij2(:,i-1,i))^(-1)*...
            (Bij2(:,i-1,i)'*((Mi(:,:,i)+Ki(:,:,i))*Dij(:,i-1,i)-(Qi(:,i)+Li(:,i)))));
        
    end
    
    Phib=[ri(:,1);euleranglei(:,1)];
    PhibZ=eye(6);
    LHS=[Ki(:,:,1)+Mi(:,:,1),PhibZ';PhibZ,zeros(6)];
    RHS=[Li(:,1)+Qi(:,1);zeros(6,1)];
    sol=LHS\RHS;
    dYi(:,1)=sol(1:6);

    for i = 1:6
%         ddqij(:,i,i+1)=-(Bij2(:,:,i,i+1)'*(Mi(:,:,i+1)+Ki(:,:,i+1))*Bij2(:,:,i,i+1))^(-1)*...
%             (Bij2(:,:,i,i+1)'*(Mi(:,:,i+1)+Ki(:,:,i+1))*Bij1(:,:,i,i+1)*dYi(:,i)+...
%             Bij2(:,:,i,i+1)'*((Mi(:,:,i+1)+Ki(:,:,i+1))*Dij(:,i,i+1)-(Qi(:,i+1)+Li(:,i+1))));
        if i == 3
            ddqij(2,i,i+1)=-(Bij2(:,i,i+1)'*(Mi(:,:,i+1)+Ki(:,:,i+1))*Bij2(:,i,i+1))^(-1)*...
                (Bij2(:,i,i+1)'*(Mi(:,:,i+1)+Ki(:,:,i+1))*Bij1(:,:,i,i+1)*dYi(:,i)+...
                Bij2(:,i,i+1)'*((Mi(:,:,i+1)+Ki(:,:,i+1))*Dij(:,i,i+1)-(Qi(:,i+1)+Li(:,i+1))));
        else
            ddqij(1,i,i+1)=-(Bij2(:,i,i+1)'*(Mi(:,:,i+1)+Ki(:,:,i+1))*Bij2(:,i,i+1))^(-1)*...
                (Bij2(:,i,i+1)'*(Mi(:,:,i+1)+Ki(:,:,i+1))*Bij1(:,:,i,i+1)*dYi(:,i)+...
                Bij2(:,i,i+1)'*((Mi(:,:,i+1)+Ki(:,:,i+1))*Dij(:,i,i+1)-(Qi(:,i+1)+Li(:,i+1))));
        end
        dqij(:,i,i+1)=dqij(:,i,i+1)+ddqij(:,i,i+1)*tstep;
        qij(:,i,i+1)=qij(:,i,i+1)+dqij(:,i,i+1)*tstep;
        
%         rij(:,i,i+1)=sij(:,i,i+1)+dij(:,i,i+1)*qij(2,i,i+1)-sji(:,i+1,i);
%         drij(:,i,i+1)=dsij(:,i,i+1)+ddij(:,i,i+1)*qij(2,i,i+1)+dij(:,i,i+1)*dqij(2,i,i+1)-dsji(:,i+1,i);
%         ri(:,i+1)=ri(:,i)+rij(:,i,i+1);
%         drijmat(:,:,i,i+1)=Vec2Mat(drij(:,i,i+1));
%         
%         dBij1(:,:,i,i+1)=[zeros(3),-drijmat(:,:,i,i+1);zeros(3),zeros(3)];
%         dBij2(:,:,i,i+1)=[-dhijmat(:,:,i,i+1)*sji(:,i+1,i)-hijmat(:,:,i,i+1)*dsji(:,i+1,i),ddij(:,i,i+1);dhij(:,i,i+1),zeros(3,1)];
%         Dij(:,i,i+1)=dBij1(:,:,i,i+1)*Yi(:,i)+dBij2(:,:,i,i+1)*dqij(:,i,i+1);
%         dYi(:,i+1)=Bij1(:,:,i,i+1)*dYi(:,i)+Bij2(:,:,i,i+1)*ddqij(:,i,i+1)+Dij(:,i,i+1);

        if i==3
            dYi(:,i+1)=Bij1(:,:,i,i+1)*dYi(:,i)+Bij2(:,i,i+1)*ddqij(2,i,i+1)+Dij(:,i,i+1);
        else
            dYi(:,i+1)=Bij1(:,:,i,i+1)*dYi(:,i)+Bij2(:,i,i+1)*ddqij(1,i,i+1)+Dij(:,i,i+1);
        end
        
    end
    euleranglei(:,1)=euleranglei(:,1)+dYi(4:6,1)*tstep;
    x7(step)=ri(1,7);
    y7(step)=ri(2,7);
end