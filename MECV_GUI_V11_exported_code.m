classdef MECV_GUI_V11_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        EnriqueCrdenasSnchezUDGUNAMLabel  matlab.ui.control.Label
        VolcanicClastCoolingRateSimulatorGUIV11Label  matlab.ui.control.Label
        TabGroup               matlab.ui.container.TabGroup
        TabFiltred             matlab.ui.container.Tab
        InterpModelListBox     matlab.ui.control.ListBox
        ModelListBoxLabel      matlab.ui.control.Label
        ImportdataButton       matlab.ui.control.Button
        SavedataSwitch         matlab.ui.control.Switch
        SavingdataSwitchLabel  matlab.ui.control.Label
        Figure1                matlab.ui.control.UIAxes
        TabSimulator           matlab.ui.container.Tab
        ImporData2             matlab.ui.control.Button
        Figure4                matlab.ui.control.UIAxes
        Figure3                matlab.ui.control.UIAxes
        Figure2                matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ImportdataButton
        function ImportdataButtonPushed(app, event)
            [file, path] = uigetfile;
             filename=[path file];
             Data=importdata(filename)
             
             Times=Data(:,1);
             Temperature=Data(:,2);

         value1 = app.InterpModelListBox.Value;
         if value1=="Linear"
                n_grade=1
          elseif value1=="Square"
              n_grade=2
          elseif value1=="Cubic"
              n_grade=3
          end 
         
 P=polyfit(Times,Temperature,n_grade);

 T_fit=polyval(P,Times);
 RR=corr(T_fit,Temperature); 
 plot(app.Figure1,Times,Temperature,Times,T_fit)
 data4saved=[Times,T_fit]
 
 value2=string(app.SavedataSwitch.Value);
 fnamedata=strcat("dataout_fit_",value1,".txt")
 
  if value2 == "Off" 
     %print("no saving data")
  elseif value2 =="On"
    cd(path)      
    save(fnamedata,"data4saved",'-ascii')
  end
        end

        % Button pushed function: ImporData2
        function ImporData2Pushed(app, event)
              [file, path] = uigetfile;
             filename=[path file];
             Data=importdata(filename);
             Times=Data(:,1); Temp_data=Data(:,2);

Solution=MECV_BySA(Data)
H_model=Solution(:,1);
Lambda_model=Solution(:,2);
Error_model=Solution(:,3);
Temp_model=Solution(:,4);


H_mle=mle(H_model)
Lambda_mle=mle(Lambda_model)


plot(app.Figure2,Times,Temp_data,Times,Temp_model)



x_index=1:12;

Dist=poisspdf(x_index,Lambda_mle(1));
plot(app.Figure3,x_index,Dist)
nbin=5;
hist(app.Figure4,H_model,nbin)

            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
function DATA=MECV_BySA(Tt_input)

 % <MECV by SIMULATED ANNEALING CODE BY HEAT TRANSFER COEFFICIENTS FROM HEAT EQUATION>
  %  Copyright (C) <2022>  <Enrique Cárdenas Sánchez>

  %  This program is free software: you can redistribute it and/or modify
  %  it under the terms of the GNU General Public License as published by
  %  the Free Software Foundation, either version 3 of the License, or
  %  (at your option) any later version.

   % This program is distributed in the hope that it will be useful,
   % but WITHOUT ANY WARRANTY; without even the implied warranty of
   % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   % GNU General Public License for more details.

   % You should have received a copy of the GNU General Public License
   % along with this program.  If not, see <https://www.gnu.org/licenses/>.


%Temperature data
T_fit=Tt_input(:,2);
% time 
t=Tt_input(:,1);

Lmax=length(T_fit);
Nd=Lmax;
index=round(linspace(1,Lmax,Nd));
tin=t(index);
Tobs=T_fit(index);
NL=7;% No. de tamaño fragmentos
%cd('D:\Respaldo de datos Camara-Termica\ProyectoTermicaColima2005')
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  lambda=seed(1), H=seed(2) Seed para comenzar la simulación.
seed=[3 20];
kont=1;
    for ii=1:Nd
    % Valor de tolerancia para la simulacón no debe de exceder de un porcentaje predeterminado..
      tol=0.005*Tobs(ii);
      % El algoritmo de Modelo de Enfriameinto de Clastos Volcanicos y
      % Simulated Annealin en un codigo llamado SimulatedAnnealin_MECV(tiempo_inicial,Temperatura_inicial,tolerancia,seed)
      [H_Best,Lambda_Best, Error, Tmodel_Best]=SimulatedAnealing_MECV(tin(ii),Tobs(ii),tol,seed);
      DE=abs(Tmodel_Best-Tobs(ii));
      %se considera el 2.5 % de tolerancia del error
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      while DE>tol
      [H_Best,Lambda_Best, Error, Tmodel_Best]=SimulatedAnealing_MECV(tin(ii),Tobs(ii),tol,seed);
      DE=abs(Tmodel_Best-Tobs(ii));
      kont=1+kont;
      end
      seed=[Lambda_Best H_Best];
      Saving_data_sol(ii,:)=[H_Best,Lambda_Best, Error, Tmodel_Best];
      DATA(ii,:)=Saving_data_sol(ii,:);
    end
end


 % <SIMULATED ANNEALING CODE BY HEAT TRANSFER COEFFICIENTS FROM HEAT EQUATION>
  %  Copyright (C) <2022>  <Enrique Cárdenas Sánchez>

  %  This program is free software: you can redistribute it and/or modify
  %  it under the terms of the GNU General Public License as published by
  %  the Free Software Foundation, either version 3 of the License, or
  %  (at your option) any later version.

   % This program is distributed in the hope that it will be useful,
   % but WITHOUT ANY WARRANTY; without even the implied warranty of
   % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   % GNU General Public License for more details.

   % You should have received a copy of the GNU General Public License
   % along with this program.  If not, see <https://www.gnu.org/licenses/>.
    

function [H_Best,Lambda_Best, Error, Tmodel_Best]=SimulatedAnealing_MECV(t,Tdata,tol,Seed)
%% data for simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     Lambda=Seed(1);
     H=Seed(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Tolerancia .... del 3.5% 
   %  tol=0.025*Tdata;

 %Definiendo variables para H
hd=200;
hc=0.026;
%Definendo valores para Lambda
ld=10; %m
lc=0.01;   %m

%%%%% Variables del simulhco recocido.
TPmax=6000;
Tp=TPmax; 
NL=7;
a_size=linspace(0.01,10,NL); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cp =@(T) 550.106+0.686*T+(-262.558)./T.^2;
K =@(T) 1./(0.3666+T*2*10^-4);
% Función de costo  
CostFunction=@(T) abs(T-Tdata);
% Funcion para h   
h=@(T,H0) H0/K(T);
  


   weight=poisspdf(1:NL,Lambda);
   VecInput=[Tdata, t, H]; 
   
 for ii=1:NL
   alpha(ii,:)=root_alpha(a_size(ii),h(Tdata,H),10);
   Temperature(ii)=Temperature_model2D(VecInput,a_size(ii),alpha(ii,:));
 end
 
 Tmodel=sum(weight.*Temperature);
 cur_cost=CostFunction(Tmodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Para esta simulhcion se considera que la maxima temperatura de experimento es constante e igual para toda las
%%% muestras. Tdata es el valor de la temperatura en el tiempo t que es menor que Tmax siempre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
                    %Finding new variables
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   new_Lambda=(ld-lc)*rand+lc;
   new_H=(hd-hc)*rand+hc;
   weight=poisspdf(1:NL,new_Lambda);
   VecInput=[Tdata, t, new_H];

 for ii=1:NL
   alpha(ii,:) = root_alpha(a_size(ii),h(Tdata,H),10);
   Temperature(ii) = Temperature_model2D(VecInput,a_size(ii),alpha(ii,:));
 end
 
 Tmodel=sum(weight.*Temperature);
 new_cost=CostFunction(Tmodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   DE=new_cost-cur_cost; % si DE es negativo curcost es mejor que newcost

if DE<0
    Lambda_Best=Lambda;
    H_Best=H;
    Error=new_cost;  
else
  new_corcost=cur_cost;
  Lambda_Best=Lambda;
  H_Best=H;
  Error=cur_cost;
end

%nuevo contador
kk=1;


%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 while Error>tol %El valor 1 representa 3.5% de grhco de erro de la Chicuhcrhca


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  choosing new solution %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Lambda=(ld-lc)*rand+lc;
   H=(hd-hc)*rand+hc;  
   weight=poisspdf(1:NL,Lambda);
   VecInput=[Tdata, t, H];

   for ii=1:NL
       alpha(ii,:)=root_alpha(a_size(ii),h(Tdata,H),10);
       Temperature(ii)=Temperature_model2D(VecInput,a_size(ii),alpha(ii,:));
   end
 
  Tmodel=sum(weight.*Temperature);
  cur_cost=CostFunction(Tmodel);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
     
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if (DE<0) 
      Lambda_Best=Lambda;
      H_Best=H; 
      Error=new_cost;  
   end
      Fa=(1-exp(-DE/TPmax))/2;
     % Existe una probabilidad de aceptar la última solución.
   if (Fa>rand/TPmax)
      Lambda_Best=Lambda;
      H_Best=H; 
      Error=cur_cost; 
      %kk=kk+1; 
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      Tp=Tp*0.98;%1/log(1+ii);
   if new_cost<cur_cost
      Lambda_Best=Lambda;
      H_Best=H;
      Error=new_cost;  
   else
      new_cost=cur_cost;
      Lambda_Best=Lambda;
      H_Best=H;
      Error=cur_cost;
   end
  kk=kk+1;
 end
  Tmodel_Best=Tmodel;       
end




 function temperature = Temperature_model2D(VecInput,a,alpha)

% Labeling variables
Tin=VecInput(1); % Kelvin Initial temperature
t=VecInput(2); % Secons
H=VecInput(3); % Value
r=0.095*a; % m Size of clast

Cp = 550.106+0.686*Tin+(-262.558)./Tin.^2;
K = 1./(0.3666+Tin*2*10^-4);

   rho=2061.9; % Valor promedio de las muestras analizadas
   kapp=K/(rho*Cp);   
   h=H/K;
  summ=0;
   
   % variables
   % Tin temperatura inicial
   % a= tamaño de muestra
   % r= distancia radial del termopar
   % alpha raices de la ecuacion transcedental
   
 for jj=1:length(alpha)
	Coeff=(a.^2.*alpha(jj).^2+(a.*h-1).^2)./(alpha(jj).^2.*(a.^2.*alpha(jj).^2+h.*a.*(a.*h-1)));
    Temp=exp(-kapp.*t.*alpha(jj).^2).*Coeff.*sin(r.*alpha(jj)).*sin(a.*alpha(jj));  
   summ=summ+Temp;
 end  
 temperature=2.*h.*Tin.*(1./r).*summ;
 end

% Encuentra las n raices para a y h
 function alpha = root_alpha (a,h,n)
aa=0.95*pi/a;
alpha=linspace(1,n);
funy =@(x)  a*x*cot(a*x)+a*h-1;
  kk=1;
     for N=1:n
         x0=N*aa;
        [x, ~] = fzero(funy,x0);
        alpha(kk)=x; 
%        ffval(kk)=fval;
        kk=kk+1;
     end   
 end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 640 399];

            % Create TabFiltred
            app.TabFiltred = uitab(app.TabGroup);
            app.TabFiltred.Title = 'Filtred Data';

            % Create Figure1
            app.Figure1 = uiaxes(app.TabFiltred);
            title(app.Figure1, 'Title')
            xlabel(app.Figure1, 'X')
            ylabel(app.Figure1, 'Y')
            zlabel(app.Figure1, 'Z')
            app.Figure1.Position = [318 149 300 185];

            % Create SavingdataSwitchLabel
            app.SavingdataSwitchLabel = uilabel(app.TabFiltred);
            app.SavingdataSwitchLabel.HorizontalAlignment = 'center';
            app.SavingdataSwitchLabel.Position = [121 186 69 22];
            app.SavingdataSwitchLabel.Text = 'Saving data';

            % Create SavedataSwitch
            app.SavedataSwitch = uiswitch(app.TabFiltred, 'slider');
            app.SavedataSwitch.Position = [132 223 45 20];
            app.SavedataSwitch.Value = 'On';

            % Create ImportdataButton
            app.ImportdataButton = uibutton(app.TabFiltred, 'push');
            app.ImportdataButton.ButtonPushedFcn = createCallbackFcn(app, @ImportdataButtonPushed, true);
            app.ImportdataButton.Icon = 'importing.png';
            app.ImportdataButton.Position = [105 128 100 22];
            app.ImportdataButton.Text = 'Import data';

            % Create ModelListBoxLabel
            app.ModelListBoxLabel = uilabel(app.TabFiltred);
            app.ModelListBoxLabel.HorizontalAlignment = 'right';
            app.ModelListBoxLabel.Position = [71 328 38 22];
            app.ModelListBoxLabel.Text = 'Model';

            % Create InterpModelListBox
            app.InterpModelListBox = uilistbox(app.TabFiltred);
            app.InterpModelListBox.Items = {'Linear', 'Square', 'Cubic', ''};
            app.InterpModelListBox.Position = [124 278 100 74];
            app.InterpModelListBox.Value = 'Linear';

            % Create TabSimulator
            app.TabSimulator = uitab(app.TabGroup);
            app.TabSimulator.Title = 'VCCRS ';

            % Create Figure2
            app.Figure2 = uiaxes(app.TabSimulator);
            title(app.Figure2, 'Model vs Data')
            xlabel(app.Figure2, 'Time (s)')
            ylabel(app.Figure2, 'Temperature (K)')
            zlabel(app.Figure2, 'Z')
            app.Figure2.Position = [181 197 278 166];

            % Create Figure3
            app.Figure3 = uiaxes(app.TabSimulator);
            title(app.Figure3, 'Grain-size Distribution')
            xlabel(app.Figure3, 'Mean Diameter (m)')
            ylabel(app.Figure3, 'Wt (%)')
            zlabel(app.Figure3, 'Z')
            app.Figure3.Position = [39 19 280 164];

            % Create Figure4
            app.Figure4 = uiaxes(app.TabSimulator);
            title(app.Figure4, 'Histogram')
            xlabel(app.Figure4, 'Diameter (m)')
            ylabel(app.Figure4, 'Wt (%)')
            zlabel(app.Figure4, 'Z')
            app.Figure4.Position = [338 15 280 172];

            % Create ImporData2
            app.ImporData2 = uibutton(app.TabSimulator, 'push');
            app.ImporData2.ButtonPushedFcn = createCallbackFcn(app, @ImporData2Pushed, true);
            app.ImporData2.Icon = 'importing.png';
            app.ImporData2.Position = [39 269 100 22];
            app.ImporData2.Text = 'Import data';

            % Create VolcanicClastCoolingRateSimulatorGUIV11Label
            app.VolcanicClastCoolingRateSimulatorGUIV11Label = uilabel(app.UIFigure);
            app.VolcanicClastCoolingRateSimulatorGUIV11Label.FontSize = 16;
            app.VolcanicClastCoolingRateSimulatorGUIV11Label.FontWeight = 'bold';
            app.VolcanicClastCoolingRateSimulatorGUIV11Label.Position = [112 445 383 22];
            app.VolcanicClastCoolingRateSimulatorGUIV11Label.Text = 'Volcanic Clast Cooling Rate Simulator GUI  V. 1.1';

            % Create EnriqueCrdenasSnchezUDGUNAMLabel
            app.EnriqueCrdenasSnchezUDGUNAMLabel = uilabel(app.UIFigure);
            app.EnriqueCrdenasSnchezUDGUNAMLabel.Position = [168 412 227 22];
            app.EnriqueCrdenasSnchezUDGUNAMLabel.Text = 'Enrique Cárdenas Sánchez, UDG-UNAM';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MECV_GUI_V11_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end