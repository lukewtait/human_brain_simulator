function xdot = JRode(t,x,pars,W)
  
vp = x(1:pars.N) ; 
ve = x(pars.N+1:2*pars.N) ; 
vi = x(2*pars.N+1:3*pars.N) ; 
zp = x(3*pars.N+1:4*pars.N) ; 
ze = x(4*pars.N+1:5*pars.N) ;
zi = x(5*pars.N+1:6*pars.N) ; 
  
vpdot = zp ; 
vedot = ze ; 
vidot = zi ; 
  
zpdot = pars.A.*pars.a.*Sigm(ve-vi)-2*pars.a.*zp-pars.a.^2.*vp ; 
zedot = pars.A.*pars.a.*(pars.p(t)+0.8*pars.C.*Sigm(pars.C.*vp) + W*vp)-2*pars.a.*ze-pars.a.^2.*ve ; 
zidot = pars.B.*pars.b.*(0.25*pars.C.*Sigm(0.25*pars.C.*vp))-2*pars.b.*zi-pars.b.^2.*vi ; 
  
xdot = [vpdot ; vedot ; vidot ; zpdot ; zedot ; zidot] ; 
  
function S = Sigm(v)
S = 2*pars.e0./(1+exp(pars.r*(pars.v0-v))) ; 
end
  
end