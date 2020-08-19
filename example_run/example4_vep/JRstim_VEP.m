function p = JRstim_VEP(t,tdelay)

if(t<=tdelay)
    st=0 ;
else
    st=0.5*(((t-tdelay)/0.005)^7)*exp(-(t-tdelay)/0.005) ;
end
    
p = 220 + [st;0;0] ;