(a+b)^2;
(a+b)^3;
(a+b)^4;
(a+b)^5;
(a+b)^6;
c:=(a+b)^6;
d:=(a+b)^3;
c/d;
a:=sin (x);
b:=cos (x);
c:=a^2+b^2;
d:=a/b;
f:=tan (x)-d;
% Trigo Rules;
trig1:={sin(~x)^2=>(1-cos(x)^2)};
let trig1;
trig2:={tan (~x)=>(sin (x)/cos (x))};
let trig2;
% now with rules;
c;
d;
f;
end;