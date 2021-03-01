% calculation of metric tensor (4-dim space-time)

off echo;
on revpri;
n:=4;

array f(n+1), dfdt(n+1), dfdl0(n+1), dfdl1(n+1), dfdl2(n+1)$

array g(n,n)$

operator x$
x(0):=t; x(1):=lambda0; x(2):=lambda1; x(3):=lambda2;

% rules
trig1:={sin(~x)^2=>(1-cos(x)^2)}$ let trig1$

% procedures
procedure scalprod(a,b); begin
integer n; n:=first(length(a))-1; 
result:=for i:=0:n-1 sum a(i)*b(i);
return result end;

procedure showmatrix(mm);begin integer m,n;l:=length(mm);m:=first(l)-1;n:=second(l)-1;matrix hm(m,n);for i:=0:m-1 do for j:=0:n-1 do hm(i+1,j+1):=mm(i,j); write hm end;

procedure showvector(vv);begin integer n;n:=first(length(vv))-1; matrix hv(n,1); for i:=0:n-1 do hv(i+1,1):=vv(i); write hv end;

% current radius
a:=a0*sqrt(1-t^2);

% surface of hyper sphere in t and lambda

f(0):=a*cos(lambda0)*cos(lambda1)*cos(lambda2);
f(1):=a*cos(lambda0)*cos(lambda1)*sin(lambda2);
f(2):=a*cos(lambda0)*sin(lambda1);
f(3):=a*sin(lambda0);
f(4):=a0*t;

for i:=0:n do dfdt(i):=df(f(i),x(0))$
for i:=0:n do dfdl0(i):=df(f(i),x(1))$
for i:=0:n do dfdl1(i):=df(f(i),x(2))$
for i:=0:n do dfdl2(i):=df(f(i),x(3))$

g(0,0):=scalprod(dfdt,dfdt)$
g(0,1):=scalprod(dfdt,dfdl0)$
g(0,2):=scalprod(dfdt,dfdl1)$
g(0,3):=scalprod(dfdt,dfdl2)$

g(1,0):=scalprod(dfdl0,dfdt)$
g(1,1):=scalprod(dfdl0,dfdl0)$
g(1,2):=scalprod(dfdl0,dfdl1)$
g(1,3):=scalprod(dfdl0,dfdl2)$

g(2,0):=scalprod(dfdl1,dfdt)$
g(2,1):=scalprod(dfdl1,dfdl0)$
g(2,2):=scalprod(dfdl1,dfdl1)$
g(1,3):=scalprod(dfdl1,dfdl2)$

g(3,0):=scalprod(dfdl2,dfdt)$
g(3,1):=scalprod(dfdl2,dfdl0)$
g(3,2):=scalprod(dfdl2,dfdl1)$
g(3,3):=scalprod(dfdl2,dfdl2)$

write "f = "; showvector(f);
write "df/dt = "; showvector(dfdt);
write "g = "; showmatrix(g);

off revpri;
on echo;
end;
