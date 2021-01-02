% potrebni toolboxi optimization, mapping
%% ustvarjanje podatkov
%N = 20;
%x1 = 10 * rand(N,1);
%y1 = 10 * rand(N,1);
%tocke = [x1,y1];

N = 4;
x1 = [0,10,5,5]';
y1 = [5,5,-2,12]';
tocke = [x1,y1];

%% Sestavljanje matrike, ki preverja ali se povezave sekajo
%seka = [];
%vrstica = 1;
%stolpec = 1;
%for i = 1:N
%    for j = 1:N
%        stolpec = 1;
%        for p = 1:N
%            for r = 1:N
%                presecisce = polyxpoly([tocke(i,1),tocke(i,2)],[tocke(j,1),tocke(j,2)],[tocke(p,1),tocke(p,2)],[tocke(r,1),tocke(r,2)]);
%                seka(vrstica,stolpec) = double(isempty(presecisce));
%                stolpec = stolpec + 1;
%                [vrstica stolpec];
%            end
%        end
%        vrstica = vrstica + 1;
%    end
%end
%seka
%
%sekanje = seka;
%for i = 0:N-1
%    sekanje((1 + N*i):N*(i+1),(1 + N*(i)):N*(i+1)) = ones(N,N); % onemogoci diagonalne bloke
%    sekanje(1 + (N+1)*i,:) = ones(1,N*N); % onemogoci vrstice, kjer je primerjava z r_ii
%    sekanje(1 + N*i:N*(i+1),(i+1):N:(i+1)+ N*(N-1)) = ones(N,N);
%    sekanje(:,1 + (N+1)*i) = ones(1,N*N); % onemogoci stolpce, kjer je primerjava z r_ii
%end
%for j = 1:N*N
%    i = mod(j,N);
%    if i == 0
%        i = 4;
%    end
%   sekanje(j,1 + (i-1)*N : (i)*N) = ones(1,N); % onemogoci vrstice, kjer je r_ij = r_jx
%   sekanje(j,i:N:N*(N-1)+i) = ones(1,N); % onemogoci vrstice, kjer je r_ij = r_xj
%end
%
%sekanje

%%
razdalje = zeros(N,N); % matrika razdalj med tockami
for i = 1 : N
    for j = 1 : N
        razdalje(i,j) = pdist([tocke(i,:);tocke(j,:)],'euclidean');
    end
end

%% racunanje

% ce racunamo min, damo po diagonali nekaj vecjega drugace pustimo na miru
% in zakomentiramo for loop
for i = 1:N
    razdalje(i,i) = max(max(razdalje))*10; % to je nekaj vecjega
end

razdalje2 = zeros(N*N,1); % vektor vseh razdalj, vzame razdalje iz matrike "razdalje" po vrsticah NxN -> N^2x1
counter = 1;
for i = 1 : N 
    for j = 1 : N
        razdalje2(counter,1) = razdalje(i,j);
        counter = counter + 1;
    end
end

f = razdalje2%*(-1); % ce iscemo maksimum pomnozimo z (-1); to so koeficienti s katerimi mnozimo y_ij
intcon = size(razdalje2,1); % nastavimo stevilo spremenljivk, ki jih iscemo
Aeq = zeros(N ^ 2,N ^ 2); % matrika kjer bom vpisal enacbe, kjer imamo enakost
beq = zeros(N ^ 2,1); % Aeq * X = Beq
lb = zeros(size(razdalje2,1),1); % spodnja meja spremenljivk je 0
ub = ones(size(razdalje2,1),1); % zgornja meja spremelnjivk je 1

% enacbe,da gre iz vsakega vozla le ena povezava
for i = 1 : N
    Aeq(i,N * (i - 1) + 1 : N * i) = 1;
    beq(i,1) = 1;
end
% enacbe, da gre v vsak vozel samo ena povezava
for i = 1 : N
   Aeq(N + i,i:N:N*(N-1)+i) = 1;
   beq(N + i,1) = 1;
end

% enacbe da se tocke povezujejo v pare
counter = 1; % N * (N-1) / 2 enacb
for j = 2:N
    for d = 0 : N-j
        Aeq(2 * N + counter,j + (N + 1) * d) = 1;
        Aeq(2 * N + counter,j + (N + 1) * d + (N-1) * (j-1)) = -1;
        beq(2 * N + counter,1) = 0;
        counter = counter + 1;
    end
    
end

%dodajanje pogoja, da se povezave ne sekajo (ne dela uredu)
%c = N*(3+N)/2 + 1;
%Aeq(c:c + N^2 - 1,:) = sekanje;
%beq(c:c + N^2 - 1,:) = 

%% resitev CLP

options = optimoptions('intlinprog','IntegerTolerance',1e-06,'RelativeGapTolerance',1e-15,'LPOptimalityTolerance',1e-10,'CutMaxIterations',50,'IntegerPreprocess','advanced')
resitev = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options)

%% graficni prikaz
matrikca = reshape(resitev,N,N);
for i = 1:N
    for j = i+1 : N
        if matrikca(i,j) == 1
            plot(tocke([i,j],1),tocke([i,j],2))
            hold on
        end
    end
end
plot(x1,y1,'o')
