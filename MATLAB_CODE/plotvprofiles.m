 function vcross = plotvprofiles(vname,timesteps)
figure(8);

clear vcross;
for i=1:length(timesteps)
    v0 = mgetfieldmpi3d([vname '0'],timesteps(i));
    v1 = mgetfieldmpi3d([vname '1'],timesteps(i));
    v2 = mgetfieldmpi3d([vname '2'],timesteps(i));
    vmag = sqrt(v0.^2 + v1.^2 + v2.^2);
    vcross(:,i)=mean(vmag)';
end
surf(vcross)