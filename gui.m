function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 29-May-2014 20:16:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in input_sig.
function input_sig_Callback(hObject, eventdata, handles)
% hObject    handle to input_sig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'stream',RandBitStream(100));
stem(getappdata(0,'stream'));

% --- Executes on button press in const_d.
function const_d_Callback(hObject, eventdata, handles)
% hObject    handle to const_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'M',16);
[ map, complex_constell ]=gray(getappdata(0,'M')); %constellation
setappdata(0,'mapped',map2gray(getappdata(0,'symbols'),map));
fig=plot(real(complex_constell),imag(complex_constell),'b.');
xlabel('In-Phase');ylabel('Quadrature');
grid on;
hold on;



% --- Executes on button press in symbols.
function symbols_Callback(hObject, eventdata, handles)
% hObject    handle to symbols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ map, complex_constell ]=gray(getappdata(0,'M'));
setappdata(0,'grouped_bits',bitgrp(getappdata(0,'stream'), getappdata(0,'M')));
setappdata(0,'symbols',binary2dec(getappdata(0,'grouped_bits')));
mapped=map2gray(getappdata(0,'symbols'), map);
fig = scatterplot(complex_constell);
scatterplot(mapped(:,2:3), 1, 0,'ro',fig)


% --- Executes on button press in noise.
function noise_Callback(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ map, complex_constell ]=gray(getappdata(0,'M'));
mapped=map2gray(getappdata(0,'symbols'), map);

setappdata(0,'k',log2(getappdata(0,'M')));

EbNo = 10;
snr = EbNo + 10*log10(getappdata(0,'k') - 10*log10(1)); 
fig = scatterplot(complex_constell);
noise=randn(size(mapped(:,2:3))); % random noise generation
constant=std(mapped(:,2:3))/(std(noise)*10^(snr/20));
setappdata(0,'receivedSignal',mapped(:,2:3) + noise*constant); %output of transmitter
setappdata(0,'noise1',noise*constant);
scatterplot(getappdata(0,'receivedSignal'), 1, 0, 'g.', fig);



% --- Executes on button press in recovery.
function recovery_Callback(hObject, eventdata, handles)
% hObject    handle to recovery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
recovered_constellation=rec_constell(getappdata(0,'receivedSignal'), map);

scatterplot((2.*round((recovered_constellation+1)/2)-1), 1, 0, 'ko', fig);

recovered_symbols=gray2symbols(recovered_constellation, map);
recovered_bits=symbols2bits(recovered_symbols);

[bit_errors, ber]=ber_calc(recovered_bits, getappdata(0,'stream'));