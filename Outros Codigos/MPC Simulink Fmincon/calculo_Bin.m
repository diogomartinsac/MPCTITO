function Bin = calculo_Bin(uk_1,umax,umin,Hc,Itil)


Bin =[repmat(umax,Hc,1) - Itil*uk_1; Itil*uk_1 - repmat(umin,Hc,1)];

end
