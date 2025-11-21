function animate(t,x,p)

q1 = x(1);
q2 = x(2);

length = p.l1+p.l2;
axis equal; axis(length*[-1,1,-1,1]); axis off;

p1 = [cos(q1); sin(q1)]*p.l1;
p2 = p1 + [cos(q2+q1); sin(q2+q1)]*p.l2;

pos = [[0;0],p1,p2];

plot(0,0,'ks','MarkerSize',25,'LineWidth',4)
plot(pos(1,:),pos(2,:),'Color',[0.1, 0.8, 0.1],'LineWidth',4)
plot(pos(1,:),pos(2,:),'k.','MarkerSize',50)

title(sprintf('Penjubot Animation,  t = %6.4f', t));



end