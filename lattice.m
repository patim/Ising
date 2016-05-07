classdef lattice < handle
    properties
        m; % lattice dimension
        data; % spins
        J; % coupling between neighboring cells
        kb; % Boltzman constant
        T; % temperature
    end
    methods
        function obj = lattice(m, T)
            % Constructor function sets up the lattice
            % gives it initial values
            % computes the energy and magnetization
            % and initializes any variables needed
            % to collect the data statistics
            obj.m = m;
            obj.data = zeros(m,m);
            obj.J = 1;
            obj.kb = 1;
            obj.T = T;
            
            for i=1:m
                for j=1:m
                    if(obj.coin()==1)
                        obj.data(i, j) = 1; 
                    else
                        obj.data(i, j) =-1;
                    end
                end
            end
            
%             obj.data(1, 1) = -1;
%             obj.data(1, 2) = 1;
%             obj.data(2, 2) = -1;
%             obj.data(2, 1) = 1;            
        end

        function mag = magnetization(obj)
            mag = 0.0;
            m = obj.m;
            for i=1:m
                for j=1:m
                    mag = mag + obj.data(i, j);
                end
            end
        end

        %
        function e = energy(obj)
            % Computing the energy of the lattice
            e = 0.0;
            m = obj.m;
            for i=1:m
                for j=1:m
                    Sij = obj.data(i, j);
                    
                    % implementing boundary conditions
                    if(i ~= m)
                        Si1 = obj.data(i+1,j);
                    else
                        Si1 = obj.data(1, j);
                    end
                    
                    if(j ~= m)                    
                        Sj1 = obj.data(i,j+1);
                    else
                        Sj1 = obj.data(i,1);
                    end
                    
                    e = e - Sij*(Si1 + Sj1);
                end
            end
        end 
        
        % Energy per spin
        function E = Energy(obj)
            N = (obj.m)^2;
            E = obj.energy()/N;
        end
        
        % Fair coin
        function p = coin(obj)
            if rand < 0.5
                p=1;
            else
                p=0;
            end
        end
        %
        function sweep(obj,N)                       
            % Implementing the Metropolis Algorithm for N random points
            for k=1:N
                i = ceil(obj.m*rand());
                j = ceil(obj.m*rand());
                
                % flipping spin
                %e_old = obj.energy();
                
                %e_new = obj.energy();
                %dE2 = e_new - e_old;
                
                dE = obj.DeltaEnergy(i, j);                
                if(dE > 0)
                    if( exp(-dE/obj.kb/obj.T) > rand())
                        % accepting the state
                        obj.data(i, j) = - obj.data(i, j);
                    end
                else
                    % keeping the state if dE is negative
                    obj.data(i, j) = - obj.data(i, j);
                end
            end
        end
        %
        function de = DeltaEnergy(obj,i,j)
            % Computing the change in the energy  of the 
            % lattice if the spin at location (i,j) flips
            
            Sij = obj.data(i,j);
            m = obj.m;
            % implementing boundary conditions
            if(i ~= m)
                Sip1j = obj.data(i+1,j);
            else
                Sip1j = obj.data(1, j);
            end
                    
            if(j ~= m)                    
                Sijp1 = obj.data(i,j+1);
            else
                Sijp1 = obj.data(i,1);
            end
            
            if(i ~= 1)
                Sim1j = obj.data(i-1,j);
            else
                Sim1j = obj.data(m, j);
            end
                    
            if(j ~= 1)                    
                Sijm1 = obj.data(i,j-1);
            else
                Sijm1 = obj.data(i,m);
            end       
            de = 2*obj.J*Sij*(Sip1j + Sim1j + Sijp1 + Sijm1);
        end
        %
        function stat(obj)
            % Accumulate the measurement statistics
        end

        % Ns - number of sweeps
        function [M,MM,E,EE] = CollectData(obj, Ns)
            % Computing and returning the averages
            N = (obj.m)^2;
            
            % performing sweeps to reach and equilibrium
            for i=1:1000
                obj.sweep(N);
            end

            Msum=0;
            MMsum=0;
            Esum=0;
            EEsum=0;
            for i=1:Ns
                obj.sweep(N);
                
                Mag = obj.magnetization();
                En = obj.energy();
                Msum = Msum + Mag;
                Esum = Esum + En;
                
                MMsum = MMsum + Mag^2;
                EEsum = EEsum + En^2;
            end
            
            % Making the avergages per spin
            M = Msum/N/Ns;
            E = Esum/N/Ns;
            MM = MMsum/N^2/Ns;
            EE = EEsum/N^2/Ns;
            
        end
        %
        function image(obj)
            % Display an image of the lattice
            % Assumes data contains +/- 1
            img = uint8(floor(0.5*(obj.data+1.0)));
            colormap([0 0 1;1 0 0]);
            image(img);
            pause(0.05);
        end
    end
end