function [] = Main_STP_plots_astro(isi, sensor)


% ---- plots for one stimulation frequency -----

% http://dgleich.github.io/hq-matlab-figs/
% 
% isi=200;
% 
% sensor='Markov';

h=figure(isi+1);
h.Position = [100 200 800 600];

Hertz=fix(1000/isi);

clf;

mainTitle = sprintf('%d Hz, %s Calcium Sensor, Astrocyte', Hertz,sensor);
sgtitle(mainTitle);

B.c    = load('csv/bc.csv');
B.Gsyn = load('csv/bg.csv');

A.aip3  = load('csv/a_ip3.csv');
A.ca    = load('csv/a_ca.csv');
A.Gsyn  = load('csv/a_Gsyn.csv');

S.Vm  = load('csv/s_Vm.csv');

tme=1:length(A.ca);
tme=tme./(1000/0.05);

fsz=8;

% length(tme)
% length(B.c)

subplot(5,1,1);
plot(tme, B.c, 'b');
t=title('Bouton Ca2+');
set(t,'FontSize', fsz);
set(gca,'FontSize', fsz);
ylabel('nM'); 
% xlabel('time (s)');

subplot(5,1,2);
plot(tme, B.Gsyn, 'g');
t=title('Synaptic [Glu]');
set(t,'FontSize', fsz);
set(gca,'FontSize', fsz);
ylabel('nM'); 
% xlabel('time (s)');

subplot(5,1,3);
plot(tme, A.aip3, 'r');
t=title('Astro [IP3]');
set(t,'FontSize', fsz);
set(gca,'FontSize', fsz);
ylabel('nM'); 
% xlabel('time (s)');

subplot(5,1,4);
plot(tme, A.ca, 'b');
t=title('Astro [Ca^{2+}]');
set(t,'FontSize', fsz);
set(gca,'FontSize', fsz);
ylabel('nM'); 
% xlabel('time (s)');

subplot(5,1,5);
plot(tme, A.Gsyn, 'g');
t=title('Extrasynaptic [Glu]');
set(t,'FontSize', fsz);
set(gca,'FontSize', fsz);
ylabel('mM'); 
xlabel('time (s)');

filename = sprintf('png/%dHz-sensor-%s-astro.png', Hertz, sensor);

saveas(h, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=figure(isi+2);
h.Position = [100 200 800 600];

plot(tme, S.Vm, 'k');

tt = sprintf('%d Hz, Spine Vm', Hertz);

t=title(tt);
set(t,'FontSize', fsz);
set(gca,'FontSize', fsz);
ylabel('mV'); 
xlabel('time (s)');

filename = sprintf('png/%dHz-sensor-%s-spine.png', Hertz, sensor);

saveas(h, filename);

