function [B, B1, B2, Bz, B1z, B2z, B1zor, Br, B1r, B2r,B1ror, Bor, B1or, B2or,...
    P, P1, P2, Por, P1or, P2or] = chebyshevBasis(m,rpoints, r)

%Matrix used to prepare even B matrix
Bd = zeros(2*m-1,rpoints); 
Bd(1,:) = ones(1,rpoints);
Bd(2,:) = r;

for i = 3:2*m-1
    Bd(i,:) = 2.*r.*Bd(i-1,:)-Bd(i-2,:);
end

% First derivative basis (2m-1 points)
Bd1 = zeros(2*m-1,rpoints);
Bd1(1,:) = zeros(1,rpoints);
Bd1(2,:) = ones(1,rpoints);

for i = 3:2*m-1
    Bd1(i,:) = 2.*Bd(i-1,:) + 2.*r.*Bd1(i-1,:)-Bd1(i-2,:);
end

% Second derivative (2m-1 points)
Bd2 = zeros(2*m-1,rpoints);
Bd2(1:2,:) = zeros(2,rpoints);

for i = 3:2*m-1
    Bd2(i,:) = 4.*Bd1(i-1,:)+2.*r.*Bd2(i-1,:)-Bd2(i-2,:);
end

% r*T(r)
Bdr = zeros(2*m-1,rpoints);
Bd1r = zeros(2*m-1, rpoints);
Bd2r = zeros(2*m-1, rpoints);

for i = 1:rpoints
    Bdr(:,i) = r(i).*Bd(:,i);
    Bd1r(:,i) = r(i).*Bd1(:,i)+Bd(:,i);
    Bd1ror(:,i) = 1/r(i).*Bd1r(:,i);
    Bd2r(:,i) = r(i).*Bd2(:,i)+2*Bd1(:,i);
end

% 1/r*T(r)
Bdor = zeros(2*m-1,rpoints);
Bd1or = zeros(2*m-1, rpoints);
Bd2or = zeros(2*m-1, rpoints);

for i = 1:rpoints
    Bdor(:,i) = 1/r(i).*Bd(:,i);
    Bd1or(:,i) = 1/r(i).*Bd1(:,i)-1/r(i).^2*Bd(:,i);
    Bd2or(:,i) = 1/r(i).*Bd2(:,i)-2/r(i)^2.*Bd1(:,i)-2/r(i).^3*Bd(:,i);
end

% Selecting only even order terms
for i = 1:m
    B(i,:) = Bd(2*i-1,:);
    B1(i,:) = Bd1(2*i-1,:);
    B2(i,:) = Bd2(2*i-1,:);
    
    Br(i,:) = Bdr(2*i-1,:);
    B1r(i,:) = Bd1r(2*i-1,:);
    B2r(i,:) = Bd2r(2*i-1,:);
    B1ror(i,:) = Bd1ror(2*i-1,:);
    
    Bor(i,:) = Bdor(2*i-1,:);
    B1or(i,:) = Bd1or(2*i-1,:);
    B2or(i,:) = Bd2or(2*i-1,:);
    
end

% % Basis with forced no-slip boundary condition
% % (1-r^2)T(r)
% for i = 1:rpoints
%     Bd2(:,i) = -2.*Bd(:,i) -4.*r(i).*Bd1(:,i) + (1-r(i).^2).*Bd2(:,i);
% end
% 
% for i = 1:rpoints
%     Bd1(:,i) = -2.*r(i).*Bd(:,i) + (1-r(i).^2).*Bd1(:,i);
% end
% 
% for i = 1:rpoints
%     Bd1zr(:,i) = 1./r(i).*Bd1(:,i);
% end
% 
% for i = 1:rpoints
%     Bd(:,i) = (1-r(i).^2).*Bd(:,i);
% end
% 
% Bz = zeros(m, rpoints);
% B1z = zeros(m, rpoints);
% B2z = zeros(m, rpoints);
% B1zor = zeros(m, rpoints);
% 
% for i = 1:m
%     Bz(i,:) = Bd(2*i-1,:);
%     B1z(i,:) = Bd1(2*i-1,:);
%     B2z(i,:) = Bd2(2*i-1,:);
%     B1zor(i,:) = Bd1zr(2*i-1,:);
% end

% Basis with forced no-slip boundary condition
% Bz_{2n}(r) = T_{2n}(r) - 1
% Bz_{0}(r) is not included in the basis so it must be removed and 2m+1
% points must be included to construct Bz basis of size m.
%Matrix used to prepare even Bz matrix
Bd = zeros(2*m+1,rpoints); 
Bd(1,:) = ones(1,rpoints);
Bd(2,:) = r;

for i = 3:2*m+1
    Bd(i,:) = 2.*r.*Bd(i-1,:)-Bd(i-2,:);
end

% First derivative basis B1z (2m+1 points)
Bd1 = zeros(2*m+1,rpoints);
Bd1(1,:) = zeros(1,rpoints);
Bd1(2,:) = ones(1,rpoints);

for i = 3:2*m+1
    Bd1(i,:) = 2.*Bd(i-1,:) + 2.*r.*Bd1(i-1,:)-Bd1(i-2,:);
end

for i = 1:rpoints
    Bd1zr(:,i) = 1./r(i).*Bd1(:,i);
end


% Second derivative B2z (2m+1 points)
Bd2 = zeros(2*m+1,rpoints);
Bd2(1:2,:) = zeros(2,rpoints);

for i = 3:2*m+1
    Bd2(i,:) = 4.*Bd1(i-1,:)+2.*r.*Bd2(i-1,:)-Bd2(i-2,:);
end

% Selecting only even order terms
for i = 1:m
    Bz(i,:) = Bd(2*i+1,:)-1;
    B1z(i,:) = Bd1(2*i+1,:);
    B1zor(i,:) = Bd1zr(2*i+1,:);
    B2z(i,:) = Bd2(2*i+1,:);
end


% Basis with forced homogeneous Neumann boundary condition
% P_{0}(r) = 1
% P_{2n}(r) = T_{2n}(r) - n^2/(n+1)^2 T_{2n+2}(r)
P = zeros(m, rpoints);
P(1,:) = ones(1, rpoints);
Bd(1,:) = ones(1,rpoints);
Bd(2,:) = r;
for i = 3:2*m+2
    Bd(i,:) = 2.*r.*Bd(i-1,:)-Bd(i-2,:);
end

for i = 2:m
    n = (2*i-1)/2;
    P(i,:) = Bd(2*n,:) - (n-0.5)^2/(1+n-0.5)^2*Bd(2*n+2,:);
end

% First derivative basis (2m-1 points)
P1 = zeros(m, rpoints);
P1(1,:) = zeros(1, rpoints);
Bd1 = zeros(2*m-1,rpoints);
Bd1(1,:) = zeros(1,rpoints);
Bd1(2,:) = ones(1,rpoints);

for i = 3:2*m+1
    Bd1(i,:) = 2.*Bd(i-1,:) + 2.*r.*Bd1(i-1,:)-Bd1(i-2,:);
end

for i = 2:m
    n = (2*i-1)/2;
    P1(i,:) = Bd1(2*n,:) - (n-0.5)^2/(1+n-0.5)^2*Bd1(2*n+2,:);
end

% Second derivative (2m-1 points)
P2 = zeros(m, rpoints);
P2(1,:) = zeros(1, rpoints);
Bd2 = zeros(2*m-1,rpoints);
Bd2(1:2,:) = zeros(2,rpoints);

for i = 3:2*m+1
    Bd2(i,:) = 4.*Bd1(i-1,:)+2.*r.*Bd2(i-1,:)-Bd2(i-2,:);
end

for i = 2:m
    n = (2*i-1)/2;
    P2(i,:) = Bd2(2*n,:) - (n-0.5)^2/(1+n-0.5)^2*Bd2(2*n+2,:);
end
Por = zeros(m, rpoints);
P1or = zeros(m, rpoints);
P2or = zeros(m, rpoints);

for i = 1:rpoints
    Por(:,i) = 1/r(i).*P(:,i);
    P1or(:,i) = 1/r(i).*P1(:,i)-1/r(i).^2*P(:,i);
    P2or(:,i) = 1/r(i).*P2(:,i)-2/r(i)^2.*P1(:,i)-2/r(i).^3*P(:,i);
end
