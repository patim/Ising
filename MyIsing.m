function[]= MyIsing()
T = 1.1; % temperature
l = lattice(64, T); % lattice with 65 cells

subplot(2,2,1)
l.image(); % initial plot of random distribution of spins
title(['Intial distribution of spins for T=', num2str(T)])

Ns=100;
[M, MM, E, EE]=l.CollectData(Ns); % collect data after 100 sweeps

Energy = l.Energy() % Energy per spin
Cv = (EE - E^2)/T^2 % Specific heat
chi = (MM - M^2)/T % Magnetic susceptibility
subplot(2,2,2)
l.image()
title(['Distribution of spins for T=',num2str(T),' after ', num2str(Ns), ' sweeps'])


T = 2.5; % temperature
l2 = lattice(64, T); % lattice with 65 cells
subplot(2,2,3)
l2.image(); % initial plot of random distribution of spins
title(['Intial distribution of spins for T=', num2str(T)])

[M, MM, E, EE]=l2.CollectData(Ns); % collect data after 100 sweeps

Energy2 = l2.Energy() % Energy per spin
Cv2 = (EE - E^2)/T^2 % Specific heat
chi2 = (MM - M^2)/T % Magnetic susceptibility

subplot(2,2,4)
l2.image()
title(['Distribution of spins for T=',num2str(T),' after ', num2str(Ns), ' sweeps'])

end