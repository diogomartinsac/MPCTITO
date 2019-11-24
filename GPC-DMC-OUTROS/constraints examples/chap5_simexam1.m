sim('chap5_exam1')
figure(1); clf reset
plot(tout,yout,'linewidth',2);
legend('Output','Input','Unsaturated Input');

sim('chap5_exam1b')
figure(2); clf reset
plot(tout,yout,'linewidth',2);
legend('Output','Input','Unsaturated Input');