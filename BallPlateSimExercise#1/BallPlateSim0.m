
%*************************************************************************
% BallPlateSim0.m
%
% BallPlate (x-AXIS ONLY) simulation base program for xPos- & xDot-
% estimation via a 2nd- or 3rd-ORDER MODEL & OBSERVER
%
% 2022/11/ß1 mys
%
% NOTES
%  Logfile ball position (XXXpos) measurements/estimates are scaled in [pxl]  
%  Logfile ball speed    (XXXdot) values/estimates       are scaled in [pxl/s]
%  They can be re-scaled into [m] and [m/s] via constant CM2PXL (see code).
%
% Format of logfile record:
%  1st column:   Abs. log time stamp          t_abs        in [s]
%  9th column:   Measured ball position       xPosMeas     in [pxl]
% 10th column:   Measured ball position       yPosMeas     in [pxl]
% 11th column:   Reference ball position      xPosRef      in [pxl]
% 12th column:   Reference ball position      yPosRef      in [pxl]
% 13th column:   Ball speed via num. dif.     xDotDifLog   in [pxl/s]
% 14th column:   Ball speed via num. dif.     yDotDifLog   in [pxl/s]
% 15th column:   Platform angle x             xAlfaPlatLog in [o]
% 16th column:   Platform angle y             yAlfaPlatLog in [o]
% 17th column:   Servo arm angle x            xThetServLog in [o]
% 18th column:   Servo arm angle y            yThetServLog in [o]
% 19th column:   Servo angle x PWM command    xThetPWMcmd  in [us]
% 20th column:   Servo angle y PWM command    yThetPWMcmd  in [us]
% 21st column:   Estimated ball x position    xPosEstLog   in [pxl]
% 22nd column:   Estimated ball x speed       xDotEstLog   in [pxl/s]
% Extended record option:
% +23rd column:  Estimated ball y position    yPosEstLog   in [pxl]
% +24th column:  Estimated ball y speed       yDotEstLog   in [pxl/s]
%
%   MatLab decimal number format uses '.' vs orig logfile uses ','  !! 
%   Original *.csv was converted to MatLab readable decimal format.
%
% Variable naming & suffix convention:
%  -Meas   real measured value
%  -Ref    real controller reference/setpoint value
%  -Dif    num. differentiation via backwards difference
%  -Pos    Position
%  -Dot    Position derivative == speed
%  -Est    Estimate (computed via Observer or KF)
%  -Obs    Estimate computed via Observer
%  -KF     Estimate computed via Kalman Filter
%
%  -Log    from log file == recorded from real system
%  -Sim    computed via this MATLAB simulation, to be compared with -Log values
%
% Changes
% 161025 Initial version, read 'BallPlateLog.cvs' & plot selected columns
% 171101 Detailed inline documentation of logfile record, added comments
% 181129 Line 206 changed sign of xThetServSim_v = - THETA_SERVO_MAX * ...
% 191025 Line 111 set PWM_X_OFFSET = 0 to compensate calibration offset
% 221101 Legends added to plots
%*************************************************************************

clear all;
close all;
clc;

% Add subfolders to MATLAB path so script can access data and helper functions
addpath('src');
addpath('data');

%*************************************************************************
% PROGRAM CONFIGURATION & SIMULATION CONTROL SECTION
%*************************************************************************

LogFile   = 'BallPlateLogData.csv';
kMax      =  300;         % # of records to be read from *.csv logfile
%kMax      =  500;         % # of records to be read from *.csv logfile

% Common/general simulation parameters
T         = 50.0e-3;      % DEFAULT sample period @ 20Hz frame rate = 50ms
kSample_v = 0:kMax-1;     % simulation sample index vector
SimTime_v = T*kSample_v;  % simulation time vector

sObsPole0 = -10.0;       % observer pole assignment in S-domain
sObsPole1 = -20.0;
%sObsPole2 = -50.0;       % only used in 3rd order system


% Various flags/switches to control program behaviour 
bLogFileFormat24Columns = 1; % if 1: extended logfile record +yEst +yDotEst in cols 23&24
bUseTrueSineTheta       = 0; % if 1: compute num. precise AlfaPlatSim with sin(ThetaServo)
bPlotBasicLogFileData   = 1; % if 1: plot real system / basic logfile data
bPlotSimulatedVsLogData = 1; % if 1: plot simulation vs real system data
bPlotOpenLoopStepTest = 1; % if 1: plot Part I open-loop step response


%*************************************************************************
% INITIALIZATION SECTION
%*************************************************************************
% Elementary numeric values / constants / scaling factors

T22       = T*T/2;        % term used in time discrete model equations

CSQRT2    = sqrt(2.0);
CSQRT3    = sqrt(3.0);
C2PI      = 2.0 * pi;
CRAD2DEG  = 180.0 / pi;   % converts [rad] -> [o]
CDEG2RAD  = pi / 180.0;   % converts [o] -> [rad]

CGEARTH   = 9.81;         % gravity constant in [m/s2]
CM2PXL    = 1380;         % camera & lens scaling factor in [pxl] per [m]
                          % for VRMagic camera (1350..1390) with 8.5mm lens
                          
% Plant model parameters
THETA_SERVO_MAX = 50;     % servo arm angle max amplitude        in [o]
PWM_CENTER   =  1500;     % PWM for servo neutral position       in [us]
PWM_AMPLTD   =   500;     % PWM amplitude for servo +/-max pos   in [us]
PWM_X_OFFSET =    45;     % PWM offset for plate horizontal pos. in [us]  -35
PWM_Y_OFFSET =     0;     % PWM offset for plate horizontal pos. in [us]  -50

mBall     = 0.003;        % ball mass        in [kg]   NOT USED/NEEDED HERE !
rBall     = 0.002;        % ball radius      in [m]    NOT USED/NEEDED HERE !
dServo    = 0.030;        % servo arm length in [m]
dPlate    = 0.200;        % plate arm attachment radius in [m]

kbb = -0.6*CGEARTH*dServo/dPlate; % model/accel. constant in SI units
kga = kbb * CM2PXL * CDEG2RAD;    % ZiSo's model/acceleration constant:
                                  % requires input (servo-) angle in [o]
                                  % returns output scaled in [pxl/s2] 

% Time CONTINUOUS 2nd order plant model:   xdot(t) = A2*x(t) + b2*u(t)
% ...
% ...
  A2=[0,1;0,0]
  b2=[0;kbb]
  c2=[1;0]
%det( Qb2 )               % observability check for CONTINUOUS system
Qb2= [c2';c2'*A2]
if (det( Qb2 ) == 0)
    warning( 'det(Q_b2) = 0: 2nd ORDER CONTINUOUS SYSTEM NOT OBSERVABLE !' ); % message/break
end


% Time DISCRETE 2nd order plant model:   x[k+1] = F2*x[k] + g2*u[k]
% ...
% ...  
F2=expm(A2*T)               % Numerical
g2=kbb.*[T22;T]
%det( Qb2 )               % observability check for TIME DISCRETE system
Qb2= [c2';c2'*F2]
if (det( Qb2 ) == 0)
    warning( 'det(Q_b2) = 0: 2nd ORDER DISCRETE SYSTEM NOT OBSERVABLE !' ); % message/break
end


% Method 2: Analytical verification using MATLAB Symbolic Math Toolbox
syms s t tau
F_s = inv(s*eye(2) - A2);          % (sI - A)^-1
F_t = ilaplace(F_s, s, t);          % Inverse Laplace transform
F2_ana = double(subs(F_t, t, T));  % Evaluate at t = T

g_t = int(F_t * b2, t, 0, T);      % Integrate e^(A*tau) * b2
g2_ana = double(g_t);              % Evaluate numerically

% Display check
disp('F2 Numerical:');   disp(F2);
disp('F2 Analytical:');  disp(F2_ana);
disp('g2 Numerical:');   disp(g2);
disp('g2 Analytical:');  disp(g2_ana);


% Time DISCRETE 2nd order observer design:  x[k+1] = F2*x[k] + g2*u[k]
z0=exp(sObsPole0*T)
z1=exp(sObsPole1*T)
%poles
q1=-(z1+z0); %(z-z0)(z-z1) coefficient (z^2-z(z1+z0)+z0z1)
q0=z0*z1;
q2=1;
% Method 1: Derived Analytic Formula
% observability gain
h2=[q1+2;
    (q0+q1+1)/T];

% Method 2: Polynomial / Ackerman Matrix Formula
% Luenberger State Observer Design
f=inv([c2';c2'*F2])*[0;1];
L=[eye(2)*q0+q1*F2+q2*F2^2]*f;          % h2= L verified

%*************************************************************************
% DATA INPUT / LOGFILE READ SECTION
%*************************************************************************

iLine = 2;    % line # to start reading (skip header & only LONG! lines) 
jCol0 = 0;    % column # to start reading

if ( bLogFileFormat24Columns == 1 ), % use extended 24 column data record
  jCol1 = 23; % last column index, extended (NEW) logfile has 24 data columns
else                                 % use 22 column data record
  jCol1 = 21; % OLDER ballplate logfile may have only 22 data columns
end

% Read 1st data record/line from logfile (i.e. skip header text line !)
LogData1_t = dlmread( LogFile , ';', iLine, jCol0, [iLine jCol0 iLine jCol1] );

% Read desired range of kMax data records into one table/matrix LogData1_t
for k = 1:kMax-1
  iLine = iLine + 2;    % current line#/record to be read from logfile
  LogRec1_t = dlmread( LogFile , ';', iLine, jCol0, [iLine jCol0 iLine jCol1] );
  %LogRec1_t            % test output of current file record 
  LogData1_t = [ LogData1_t
                 LogRec1_t  ]; % append new record to total table
end % for 

% Split up logfile table/matrix into individual variables -> column vectors XXX_v
AbsTime_v      = LogData1_t(:,1);  %  1st column:  absolute log time in [s]

xPosMeas_v     = LogData1_t(:,9);  %  9th column:  xPosMeas   in [pxl]
yPosMeas_v     = LogData1_t(:,10); % 10th column:  yPosMeas   in [pxl]

xPosRef_v      = LogData1_t(:,11); % 11th column:  xPosRef    in [pxl]
yPosRef_v      = LogData1_t(:,12); % 12th column:  yPosRef    in [pxl]

xDotDifLog_v   = LogData1_t(:,13); % 13th column:  xDotDifLog in [pxl/s]
yDotDifLog_v   = LogData1_t(:,14); % 14th column:  yDotDifLog in [pxl/s]

xAlfaPlatLog_v = LogData1_t(:,15); % 15th column:  platform angle x  in [o]
yAlfaPlatLog_v = LogData1_t(:,16); % 16th column:  platform angle y  in [o]

xThetServLog_v = LogData1_t(:,17); % 17th column:  servo arm angle x in [o]
yThetServLog_v = LogData1_t(:,18); % 18th column:  servo arm angle y in [o]

xThetPWMcmd_v = LogData1_t(:,19)-PWM_CENTER+PWM_X_OFFSET; % 19th col: xThetPWMcmd in [us]
yThetPWMcmd_v = LogData1_t(:,20)-PWM_CENTER+PWM_Y_OFFSET; % 20th col: yThetPWMcmd in [us]

xPosEstLog_v = LogData1_t(:,21);   % 21st column:  xPosEstLog    in [pxl]
xDotEstLog_v = LogData1_t(:,22);   % 22nd column:  xDotEstLog    in [pxl/s]

if ( bLogFileFormat24Columns == 1 ),
  yPosEstLog_v = LogData1_t(:,23); % +23rd column: +yPosEstLog   in [pxl]
end

if ( bLogFileFormat24Columns == 1 ),
  yDotEstLog_v = LogData1_t(:,24); % +24th column: +yDotEstLog   in [pxl/s]
end


% compute alfa platform simulation values (+/- numerically precise/true sin(Theta) )
xThetServSim_v = - THETA_SERVO_MAX * xThetPWMcmd_v / PWM_AMPLTD;      % in [o]

if ( bUseTrueSineTheta == 1 ) % use true sin(ThetaServo)
  xAlfaPlatSim_v = sin(xThetServSim_v*CDEG2RAD) * dServo / dPlate;
else                          % use linearization for sin(ThetaServo)
  xAlfaPlatSim_v = xThetServSim_v * dServo / dPlate;            % in [o]
end


if ( 0 ),
% Plot xAlfaPlatLog vs xAlfaPlatSim
  figure(2); clf;
  plot(SimTime_v, [ xThetServLog_v  xThetServSim_v ] );
  grid on;
  xlabel('SimTime [s]')
  ylabel('Theta  [o] [o]')
  legend('ThetServLog', 'ThetServSim');
  title('ThetaServLog  ThetaServSim');

  figure(4); clf;
  plot(SimTime_v, [ xAlfaPlatLog_v  xThetServSim_v  xAlfaPlatSim_v ] );
  grid on;
  xlabel('SimTime [s]')
  ylabel('Alfa/Theta  [o] [o] [o]')
  legend('AlfaPlatLog', 'ThetServSim', 'AlfaPlatSim');
  title('AlfaPlatLog  ThetaServSim AlfaPlatSim');
end

%*************************************************************************
% DATA PROCESSING & SIMULATION SECTION
%*************************************************************************

% dimension & init simulation vars / plot vectors
xDotSimLPF_v   = zeros( kMax,1 );

% 2nd order system simulation
xPosSim2Obs_v  = zeros( kMax,1 ); % xPos2 via observer simulation
xDotSim2Obs_v  = zeros( kMax,1 ); % xDot2 via observer simulation

% dimension & init 2nd order state vector before starting sim loop
xState2Obsk  = [ xPosMeas_v(1)    % initialize xObsk(1) w measured value
                 0 ];
xState2Obsk1 = [ xPosMeas_v(1)    % initialize xObsk1(1) w measured value
                 0 ];
             

% --- PART I: OPEN-LOOP STEP RESPONSE SIMULATION ---
u_step = 30 * (SimTime_v >= 1.0) ...
    - 60 * (SimTime_v >= 2.0) ...
    + 30 * (SimTime_v >= 3.0); % [deg]

xSim2_openloop = zeros(2, kMax);

for k = 1:kMax-1
    xSim2_openloop(:, k+1) = F2 * xSim2_openloop(:, k) + g2 * u_step(k);
end



% --- PART II: OBSERVER SIMULATION ---

uk     = xThetServSim_v; % scalar control input variable (here: servo arm angle Theta in [o])
y2 = xPosMeas_v;       % Measured position [pxl]
x2=zeros(kMax+1,2).';       %OR x2 = zeros(2, kMax+1);
%  Initialize observer state with initial measured position
x2(:, 1) = [xPosMeas_v(1); 0];

% MAIN SIMULATION LOOP
for k = 1:kMax            % sample index 
    % Time DISCRETE 2nd order observer design:
    %            x[k+1] = F2*x[k] + g2*u[k] + h2*( y[k] - c2'x[k] )    
    % ...if xThetServSim_v=
    x2(:,k+1)=F2*x2(:,k)+g2*uk(k)+h2*(y2(k)-(c2'*x2(:,k)));

end % MAIN SIMULATION LOOP

% 5. Store results for comparison plots
xPosSim2Obs = x2(1, 1:kMax)'; % Observer estimated position [pxl]
xDotSim2Obs = x2(2, 1:kMax)'; % Observer estimated velocity [pxl/s]

%*************************************************************************
% DATA OUTPUT / PLOTTING SECTION
%*************************************************************************


if ( 0 )  % 1/0 enables/disables the following code block
  figure(5); clf;  % Plot xPosRef/Meas vs. xThetPWMcmd
  plot(AbsTime_v, [ xPosRef_v  xPosMeas_v  xThetPWMcmd_v ] );
  grid on;
  xlabel('AbsTime [s]')
  ylabel('xPos  [pxl] [pxl]  [us]')
  title('xPosRef(blu)  xPosMeas(red)  xThetPWMcmd(yel)')
  
  figure(100); clf;
  plot(AbsTime_v, [ xPosMeas_v  yPosMeas_v ] );
 %plot(SimTime_v, [ xPosMeas_v  yPosMeas_v ] );
  grid on;
  xlabel('AbsTime [s]')
  ylabel('x/yPos  [pxl] [pxl]')
  title('xPosMeas(blu)  yPosMeas(red)')  
end

%*************************************************************************
% PART I: OPEN-LOOP STEP RESPONSE PLOT
%*************************************************************************
if ( bPlotOpenLoopStepTest == 1 )
    figure(100); clf;

    subplot(3,1,1);
    plot(SimTime_v, u_step, 'r', 'LineWidth', 1.5); 
    grid on;
    ylabel('xThetServo [^o]');
    title('Part I: Open-Loop Step Response');

    subplot(3,1,2);
    plot(SimTime_v, xSim2_openloop(1,:), 'b', 'LineWidth', 1.5); 
    grid on;
    ylabel('xPosSim [pxl]');

    subplot(3,1,3);
    plot(SimTime_v, xSim2_openloop(2,:), 'k', 'LineWidth', 1.5); 
    grid on;
    xlabel('SimTime [s]'); 
    ylabel('xDotSim [pxl/s]');
end


if ( bPlotBasicLogFileData == 1 ) % plot basic real system/logfile data
   
  % Plot REAL SYSTEM's BALL POSITION ... xPosRef/Meas vs. xPosEstLog
  figure(10); clf;
  plot(SimTime_v, [ xPosRef_v  xPosMeas_v  xPosEstLog_v ] );
  grid on;
  xlabel('SimTime [s]')
  ylabel('xPos  [pxl] [pxl] [pxl] ')
  legend('xPosRef [pxl]', 'xPosMeas [pxl]', 'xPosEst [pxl]');
  title('REAL SYSTEM  xPosRef xPosMeas  vs  xPosEst') 

  % Plot LOGDATA's BALL SPEED ... xDotDifLog/Est vs. xThetPWMcmd
  figure(20); clf;
  plot(SimTime_v, [ xThetServLog_v*10  xDotDifLog_v  xDotEstLog_v ] );
  grid on;
  xlabel('SimTime [s]')
  ylabel('xDot  [pxl/s] [pxl/s]   xThetServ x 10[o]')
  legend('xThetServo [o]', 'xDotDiff [pxl/s]', 'xDotEst [pxl/s]');
  title('REAL SYSTEM  xThetServLog xDotDiff  vs  xDotEst')
  
  figure(30); clf;
  plot(SimTime_v, [ xThetPWMcmd_v  xDotDifLog_v ] );
  grid on;
  xlabel('SimTime [s]')
  ylabel('xDot  [pxl/s]   xThetPWMcmd [us] ')
  legend('ThetPWMcmd [us]', 'xDotDiff [pxl/s]');
  title('REAL SYSTEM  xThetPWMcmd  vs  xDotDiff')

  if ( 0 )
    figure(50); clf;
    plot(SimTime_v, [ xAlfaPlatLog_v*100 xDotDifLog_v  ] );
    grid on;
    xlabel('SimTime [s]')
    ylabel('xDot  [pxl/s]   xAlfaPlat x 100[o] ')
    title('REAL SYSTEM xAlfaPlat(blu) vs xDotDiff(red) ')
  end

end % plot basic real system/logfile data


%*************************************************************************
% PART II: OBSERVER SIMULATION PLOT
%*************************************************************************

if ( bPlotSimulatedVsLogData == 1 ) % plot simulation vs real system/logfile data
    
 figure(200); clf;

    % --- Position Comparison ---
    subplot(2,1,1);
    plot(SimTime_v, xPosMeas_v, 'k.', 'MarkerSize', 5); hold on;
    plot(SimTime_v, xPosEstLog_v, 'r--', 'LineWidth', 1.5);
    plot(SimTime_v, xPosSim2Obs, 'b-', 'LineWidth', 1.2);
    grid on;
    ylabel('xPos [pxl]');
    title('Part II: 2nd Order Luenberger Observer Estimates vs Logfile Data');
    legend('xPosMeas (Measured)', 'xPosEstLog (Logfile)', 'xPosSim2Obs (Our Observer)');

    % --- Velocity Comparison ---
    subplot(2,1,2);
    plot(SimTime_v, xDotDifLog_v, 'g:', 'LineWidth', 1); hold on;
    plot(SimTime_v, xDotEstLog_v, 'r--', 'LineWidth', 1.5);
    plot(SimTime_v, xDotSim2Obs, 'b-', 'LineWidth', 1.2);
    grid on;
    xlabel('SimTime [s]');
    ylabel('xDot [pxl/s]');
    legend('xDotDifLog (Num Diff)', 'xDotEstLog (Logfile)', 'xDotSim2Obs (Our Observer)');
  
end % plot simulation vs real system/logfile data



%*************************************************************************
% End of BallPlateSim0.m
%*************************************************************************

