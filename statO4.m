function [A,Fp,Fu,Fw]=statO4(n)
mat=material;
X=netz(n);
EL=elem(X); h=abs(X(2,1)-X(1,1));
% 1.Nebenbedingung (Zustands-Gleichung)
% int phi(u)d_x^2w d_x^2v +psi(u)wv dx = int fv dx (w)
% Adjungierte Gleichung
% int phi(u)d_x^2p d_x^2v +psi(u)pv dx = int v dx (p)
% Input Elemente, Knoten
[c1,c2]=size(X);
[e1,e2]=size(EL);
A=sparse(2*c1,2*c1); % Sparce-Matrizen und Vektoren
Fw=sparse(2*c1,c1); % rhs w
Fp=sparse(2*c1,2*c1); % rhs p
Fu=sparse(c1,c1); % obj u
% Elementroutine, geht ueber alle Elemente und berechnet
% Elementarmatrizen, die in der Steifigkeitsmatrix A
% eingebaut werden. Analoge Vorgehensweise bei den rechten Seiten f_w und
% f_p
for k=1:e1
Knoten=[EL(k,1) EL(k,2)];
[A_E,fw,fp,fu]=elemm(Knoten,X,h,mat);
j=1;
for i=Knoten % pro Knoten 2 FG
Knoten2(j)=(2*i-1); j=j+1;
end
KH=[1 3];
 for j=1:2
  for i=1:2
    A(Knoten2(j):Knoten2(j)+1,Knoten2(i):Knoten2(i)+1)=...
    A(Knoten2(j):Knoten2(j)+1,Knoten2(i):Knoten2(i)+1)+...
    A_E(KH(j):KH(j)+1,KH(i):KH(i)+1);
Fw(Knoten2(j):Knoten2(j)+1,k+i-1)=Fw(Knoten2(j):Knoten2(j)+1,k+i-1)+...
fw(KH(j):KH(j)+1,i);
Fp(Knoten2(j):Knoten2(j)+1,Knoten2(i):Knoten2(i)+1)=Fp(Knoten2(j):Knoten2(j)+1,Knoten2(i):Knoten2(i)+1)+...
fp(KH(j):KH(j)+1,KH(i):KH(i)+1);
Fu(k+j-1,k+i-1)=Fu(k+j-1,k+i-1)+fu(j,i);
  end 
 end 
end
end
%a: XG*(1-XG)
%b: exp(XG)
%Routinen zum Berechnen des Finite-Elemente Gradienten:
%Routinen für die Berechnung des Gradienten der Variationsungleichung:
%----------------- List of Functions --------------------------