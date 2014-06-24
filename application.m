%%application.m
%BY PRZEMYS£AW DYMITROWSKI
%   PIOTR SZKOTAK


function varargout = application(varargin)
% APPLICATION MATLAB code for application.fig
%      APPLICATION, by itself, creates a new APPLICATION or raises the existing
%      singleton*.
%
%      H = APPLICATION returns the handle to a new APPLICATION or the handle to
%      the existing singleton*.
%
%      APPLICATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APPLICATION.M with the given input arguments.
%
%      APPLICATION('Property','Value',...) creates a new APPLICATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before application_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to application_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help application

% Last Modified by GUIDE v2.5 25-Jun-2014 00:50:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @application_OpeningFcn, ...
                   'gui_OutputFcn',  @application_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before application is made visible.
function application_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to application (see VARARGIN)

% Choose default command line output for application
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes application wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = application_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function m_text_Callback(hObject, eventdata, handles)
% hObject    handle to text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text as text
%        str2double(get(hObject,'String')) returns contents of text as a double


handles.text=str2double(get(hObject,'string'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function m_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bistream_text_Callback(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text1 as text
%        str2double(get(hObject,'String')) returns contents of text1 as a double


handles.text1=get(hObject,'string');
set(handles.text1,'String',stream);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function bistream_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_text_Callback(hObject, eventdata, handles)
% hObject    handle to n_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_text as text
%        str2double(get(hObject,'String')) returns contents of n_text as a double


handles.n_text=str2double(get(hObject,'string'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function n_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gen_bits_btn.
function gen_bits_btn_Callback(hObject, eventdata, handles)
% hObject    handle to gen_bits_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global stream
global N;
% global stream;
N=handles.n_text;   
stream = RandBitStream(N);
stream = num2str(transpose(stream)); % we convert it to string
stream = stream(~isspace(stream))'; % we remove the spaces
stream = strtrim(stream);% trim
set(handles.bistream_text,'String',stream); % we set the string of bitstream to the edit box
guidata(hObject, handles);







% --- Executes on button press in mod_btn.
function mod_btn_Callback(hObject, eventdata, handles)
% hObject    handle to mod_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global M;
global N;


N = handles.n_text;
M = handles.text;

bistream_text = get(handles.bistream_text, 'String');
input_string = bistream_text;
for n=1:length(input_string)
stream(n)=bin2dec(input_string(n));
end


k=log2(M); % bits per symbol

[stream, nb]=prepare_stream(stream, k); %zero padding if needed
axes(handles.input_bit);
stem(stream);


EbNo = 20; % to calculate snr
fc=40; %frequency of the carrier
 
% 
 fs=500;     	    % sampling frequency [Hz]
% frs=N*round(fs/N); % real sampling frequency [Hz]
% fn=frs/2;           % Nyquist frequency [Hz]
% Ts=1/frs;	        % Sampling time [s]
% t=[0:Ts:(N*Tb)-Ts]';% Time vector initialization


frs=nb*round(fs/nb); % real sampling frequency [Hz]
Ts=1/frs;	        % Sampling time [s]
Tb=1/nb;		        % Bit duration time [s]
t=[0:Ts:(nb*Tb)-Ts]';% Time vector initialization



%mapping
[ map, complex_constell ]=gray(M); %creating constellation
grouped_bits=bitgrp(stream, M); %grouping bits
symbols=binary2dec(grouped_bits); % creating symbols
mapped=map2gray(symbols, map); %mapping symbols to the constellation points
%plotting constellation 
axes(handles.const_axes)
plot(real(complex_constell),imag(complex_constell),'b.');


%splitting
[I, Q]=split_stream(mapped);



% creating modulating signals sI and sQ

len=nb/k; %length of I or Q in symbols
rep=frs/len; %number of repetitions of a single bit (to obtain square wave)

x=1:1:(len+1)*(1/Ts*k);
for i=1:len
    for j=1:.1:i+1;
        sI(x(i*rep:(i+1)*rep))=I(i);
        sQ(x(i*rep:(i+1)*rep))=Q(i);
    end
end
setappdata(0,'sI',sI(end-frs:end-1));
setappdata(0,'sQ',sQ(end-frs:end-1));

% Carriers
carrier_I=cos(2*pi*fc*t);
carrier_Q=-sin(2*pi*fc*t);

% Modulating
setappdata(0,'modulated_I',getappdata(0,'sI')'.*carrier_I);
setappdata(0,'modulated_Q',getappdata(0,'sQ')'.*carrier_Q);

% Summing I and Q
output_signal=getappdata(0,'modulated_I')+getappdata(0,'modulated_Q');

%plots 
axes(handles.in_axes)
plot(t, getappdata(0,'sI'), 'r'); title ('Modulating signal I'); ylabel ('I amplitude');  xlabel ('t[s]'); 
axes(handles.in_car_axes)
plot(t, getappdata(0,'modulated_I'), 'r'); title ('Modulated I'); ylabel ('I amplitude');  xlabel ('t[s]'); 

axes(handles.quad_axes)
plot(t, getappdata(0,'sQ'), 'g'); title ('Modulating signal Q'); ylabel ('Q amplitude'); xlabel ('t[s]'); 
axes(handles.quad_car_axes)
plot(t, getappdata(0,'modulated_Q'), 'g'); title ('Modulated Q'); ylabel ('Q amplitude');  xlabel ('t[s]');
axes(handles.out_axes)
plot(t, output_signal); title ('Output Signal - sum of modulated I and Q'); ylabel ('Amplitude');  xlabel ('t[s]'); 

%transmittiing through AWGN

%calculating SNR
snr = EbNo + 10*log10(k);

noise=randn(size(output_signal)); % random noise generation
constant=std(output_signal)/(std(noise)*10^(snr/20));
setappdata(0,'received_signal',output_signal + noise*constant); %output of transmitter

axes(handles.out_nois_axes)
plot(t,getappdata(0,'received_signal')); title ('Output Signal - transmitted through AWGN channel'); ylabel ('Amplitude');  xlabel ('t[s]');



% --- Executes on button press in demod_btn.
function demod_btn_Callback(hObject, eventdata, handles)
% hObject    handle to demod_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
 global M;
 global N;
 
 M = handles.text;
 N = handles.n_text;
 k=log2(M); % bits per symbol
 
[stream, nb]=prepare_stream(getappdata(0,'stream'), k); %zero padding if needed


EbNo = 20; % to calculate snr
fc=40; %frequency of the carrier
Tb=1/nb;		        % Bit duration time [s]
 

fs=500;     	    % sampling frequency [Hz]

frs=nb*round(fs/nb); % real sampling frequency [Hz]
Ts=1/frs;	        % Sampling time [s]
Tb=1/nb;		        % Bit duration time [s]
t=[0:Ts:(nb*Tb)-Ts]';% Time vector initialization


% Carriers
carrier_I=cos(2*pi*fc*t);
carrier_Q=-sin(2*pi*fc*t);

%demodulation

I_recovered=getappdata(0,'received_signal').*carrier_I;
Q_recovered=getappdata(0,'received_signal').*carrier_Q;


%filtration

filtered_I=lpf(I_recovered, fc, frs);
filtered_Q=lpf(Q_recovered, fc, frs);


filtered_I=2*filtered_I;
filtered_Q=2*filtered_Q;

%averaging

sym_num=nb/k;
demodulated_I=mean(reshape(filtered_I, ceil(frs/sym_num), []))';
demodulated_Q=mean(reshape(filtered_Q, ceil(frs/sym_num), []))';
% received_signal=horzcat(demodulated_I, demodulated_Q)
% mappedSignal=mapped(:,2:3) %just to display and compare with receivedSignal


% %plots, constellation
% fig = scatterplot(complex_constell); pause;
% hold on;
% scatterplot(mapped(:,2:3), 1, 0,'ro',fig); pause;
% scatterplot(receivedSignal, 1, 0, 'g.', fig); pause;
% hold off;
% close all;

%inverse mapping and symbols to bitstream conversion
% recovered_constellation=rec_constell(getappdata(0,'receivedSignal'), map);
% recovered_symbols=gray2symbols(recovered_constellation, map); 
% setappdata(0,'recovered_bits',symbols2bits(recovered_symbols));

%plots in frequency most important blue and red!

axes(handles.fft_axes)
cla(gca);
hold on;
plot_fft(getappdata(0,'sI'), frs, 'b', 'Spectrum of modulating signal sI'); 
plot_fft(getappdata(0,'modulated_I'), frs, 'g', 'Spectrum of modulated signal (sI*carrier)'); 
plot_fft(I_recovered, frs, 'y', 'Spectrum of modulated signal after multiplication by carrier at the receiver'); 
plot_fft(2*filtered_I, frs, 'r', 'Spectrum of demodulated signal after LPF'); 
hold off;

%plots time
axes(handles.in_axes);
cla(gca);hold on;
plot(t, getappdata(0,'sI'), 'b'); title('Modulating signal I'); xlabel('t[s]'); ylabel('Amplitude');
plot(t, filtered_I, 'r'); title('Recovered I'); xlabel('t[s]'); ylabel('A');
hold off;

axes(handles.quad_axes);
cla(gca); hold on;
plot(t, getappdata(0,'sQ'), 'b'); title('Modulating signal Q'); xlabel('t[s]'); ylabel('Amplitude');  
plot(t, filtered_Q, 'r'); title('Recovered Q'); xlabel('t[s]'); ylabel('A'); 
hold off;

%real ber
% [bit_errors, ber]=ber_calc(getappdata(0,'recovered_bits'), getappdata(0,'stream'));

%theoretical ber, some approximations taken // wiki
x=sqrt(3*k*EbNo/(M-1));
theoretical_ber=(4/k)*(1-1/sqrt(M))*(1/2)*erfc(x/sqrt(2));

% bit_errors
% ber
% theoretical_ber
% length(recovered_bits)
