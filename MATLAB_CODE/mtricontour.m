function mtricontour(timeStep,rep)

data0 = mgetfieldmpi3d('phi0',timeStep);
data1 = mgetfieldmpi3d('phi1',timeStep);
data2 = mgetfieldmpi3d('phi2',timeStep);

phi0 = repmat(data0,[2 1 1]);
phi1 = repmat(data1,[2 1 1]);
phi2 = repmat(data2,[2 1 1]);

temp0(:,:)=phi0(1,:,:);
temp1(:,:)=phi1(1,:,:);
temp2(:,:)=phi2(1,:,:);

surf(temp0,'FaceColor',hsv2rgb([0/3 1 1]),'EdgeColor','none')
hold on
surf(temp1,'FaceColor',hsv2rgb([1/3 1 1]),'EdgeColor','none')
surf(temp2,'FaceColor',hsv2rgb([2/3 1 1]),'EdgeColor','none')
hold off

lighting phong;
camlight left;
%material([0.5 0.5 0.5]);
material dull;
axis square tight;
zlim([0 1]);
view(3);
%view([0 90]);
pbaspect([1 1 0.20])