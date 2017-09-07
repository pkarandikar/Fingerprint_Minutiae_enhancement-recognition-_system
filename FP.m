function varargout = FP(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FP_OpeningFcn, ...
                   'gui_OutputFcn',  @FP_OutputFcn, ...
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
function FP_OpeningFcn(hObject, eventdata, handles, varargin)
if (exist('FingerPrint_database.dat')==2)
            load('FingerPrint_database.dat','-mat');
            disp('               *********   Database Info   *********               ');
          w='               *********   Database Info   *********           ';
          for ii=1:fp_number
                lunghezzanome = length(namefile_vector{ii});
                aggiuntaspazi = 30 - lunghezzanome;
                if aggiuntaspazi<0
                    aggiuntaspazi = 3;
                end
                stringavuota      = '-';
                stringaaggiungere = '';
                for jj=1:aggiuntaspazi
                    stringaaggiungere = strcat(stringavuota,stringaaggiungere);
                end
               
                set(handles.edit3,'String',fp_number);
                
            end
        else
            
            w='Database is Clear';
            
        end

handles.output = hObject;


guidata(hObject, handles);

function varargout = FP_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global pathname;
global namefile;
global img;
global selezionato;

clc;
        selezionato=0;

        [namefile,pathname]=uigetfile({'*.bmp;*.tif;*.tiff;*.jpg;*.jpeg;*.gif;*.pgm','IMAGE Files (*.bmp,*.tif,*.tiff,*.jpg,*.jpeg,*.gif,*.pgm)'},'Chose GrayScale Image');
        if namefile~=0
            [img,map]=imread(strcat(pathname,namefile));
            axes(handles.axes2);
            imshow(img);
            [oimg,fimg,bwimg,eimg,enhimg] =  fft_enhance_cubs(img);
             axes(handles.axes3);
            imshow(enhimg);
            [gimg,oimg] = orientation_image_rao(img);
             axes(handles.axes4);
             view_orientation_image(oimg);
             [oimg,fimg,bwimg,eimg,enhimg] =  fft_enhance_cubs(img);
                [xc,yc]=supercore7(enhimg);
                
                 axes(handles.axes5);
                imshow(img);
                hold on;
                plot(yc,xc,'O');
                hold off;
            selezionato=1;
        else
            w='Select a grayscale image';
            
        end
        if (any(namefile~=0) &&  length(size(img))>2)
            w='Select a grayscale image';            
            selezionato=0;
        end

function axes2_CreateFcn(hObject, eventdata, handles)

function pushbutton2_Callback(hObject, eventdata, handles)

global pathname;
global namefile;
global img;
global selezionato;
global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice matricer num_disk


n_bands=4;
h_bands=20;
n_arcs=16;
h_radius=12;
h_lato=h_radius+(n_bands*h_bands*2)+16;
if mod(h_lato,2)==0
    h_lato=h_lato-1;
end
n_sectors=n_bands*n_arcs;
num_disk=8;

matrice  = zeros(h_lato); 
matricer = zeros(h_lato); 
for ii=1:(h_lato*h_lato)
    matrice(ii) = whichsector(ii,0);
    matricer(ii)= whichsector(ii,1);
end
 if (any(namefile~=0) &&  length(size(img))>2)
            w='Select a grayscale image';
           
            selezionato=0;
        end

        if selezionato==1


            immagine=double(img);

            if isa(img,'uint8')
                graylevmax=2^8-1;
            end
            if isa(img,'uint16')
                graylevmax=2^16-1;
            end
            if isa(img,'uint32')
                graylevmax=2^32-1;
            end
            fingerprint = immagine;

            N=h_lato;

            [oimg,fimg,bwimg,eimg,enhimg] =  fft_enhance_cubs(fingerprint);
            fingerprint = enhimg;
            [YofCenter,XofCenter]=supercore7(fingerprint);
            [CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
            [NormalizedPrint,vector]=sector_norm(CroppedPrint,0,0);

            for (angle=0:1:num_disk-1)
                gabor=gabor2d_sub(angle,num_disk,0);
                ComponentPrint=conv2fft(NormalizedPrint,gabor,'same');
                [disk,vector]=sector_norm(ComponentPrint,1,0);
                finger_code1{angle+1}=vector(1:n_sectors);
            end

            
            [NormalizedPrint,vector]=sector_norm(CroppedPrint,0,1);

            for (angle=0:1:num_disk-1)
                gabor=gabor2d_sub(angle,num_disk,1);
                ComponentPrint=conv2fft(NormalizedPrint,gabor,'same');
                [disk,vector]=sector_norm(ComponentPrint,1,1);
                finger_code2{angle+1}=vector(1:n_sectors);
            end
            
            if (exist('FingerPrint_database.dat')==2)
                load('FingerPrint_database.dat','-mat');
                fp_number=fp_number+1;
                data{fp_number,1}=finger_code1;
                data{fp_number,2}=finger_code2;
                namefile_vector{fp_number} = namefile;
                path_vector{fp_number}     = pathname;
                save('FingerPrint_database.dat','data','fp_number','namefile_vector','path_vector','-append');
            else
                fp_number=1;
                data{fp_number,1}=finger_code1;
                data{fp_number,2}=finger_code2;
                namefile_vector{fp_number} = namefile;
                path_vector{fp_number}     = pathname;
                save('FingerPrint_database.dat','data','fp_number','namefile_vector','path_vector');
            end

            message=strcat('Succesfully added to database. FingerCode no. ',num2str(fp_number));
            msgbox(message,'FingerCode DataBase','help');
        end
load('FingerPrint_database.dat','-mat');  
set(handles.edit3,'String',fp_number);

function pushbutton3_Callback(hObject, eventdata, handles)

global pathname;
global namefile;
global img;
global selezionato;
n_bands=4;
h_bands=20;
n_arcs=16;
h_radius=12;
h_lato=h_radius+(n_bands*h_bands*2)+16;
if mod(h_lato,2)==0
    h_lato=h_lato-1;
end
n_sectors=n_bands*n_arcs;
num_disk=8;

matrice  = zeros(h_lato); 
matricer = zeros(h_lato); 
for ii=1:(h_lato*h_lato)
    matrice(ii) = whichsector(ii,0);
    matricer(ii)= whichsector(ii,1);
end
 if (exist('FingerPrint_database.dat')==2)
            load('FingerPrint_database.dat','-mat');

            
            if (any(namefile~=0) &&  length(size(img))>2)
                w='Select a grayscale image';
            
                selezionato=0;
            end

            if selezionato==1

                message = strcat('Image selected for fingerprint matching: ',namefile);
               
           
                message = strcat('Location: ',pathname);
                disp(message);

                immagine=double(img);

                if isa(img,'uint8')
                    graylevmax=2^8-1;
                end
                if isa(img,'uint16')
                    graylevmax=2^16-1;
                end
                if isa(img,'uint32')
                    graylevmax=2^32-1;
                end
                fingerprint = immagine;

                N=h_lato;

                [oimg,fimg,bwimg,eimg,enhimg] =  fft_enhance_cubs(fingerprint);
                fingerprint = enhimg;
                [list_of_core]= supercore7_list(fingerprint);
                results_img       = zeros(size(list_of_core,1),1);
                results_dis       = zeros(size(list_of_core,1),1);
                message = strcat('Candidates for core points: ',num2str(size(list_of_core,1)));
               
                for scan_core = 1:size(list_of_core,1)
                    message = strcat('Scanning candidate # ',num2str(scan_core));
                   

                    YofCenter = list_of_core(scan_core,1);
                    XofCenter = list_of_core(scan_core,2);
                   
                    [CroppedPrint]=cropping(XofCenter,YofCenter,fingerprint);
                    [NormalizedPrint,vector]=sector_norm(CroppedPrint,0,0);

                    
                    vettore_in=zeros(num_disk*n_sectors,1);
                    for (angle=0:1:num_disk-1)
                        gabor=gabor2d_sub(angle,num_disk,0);
                        ComponentPrint=conv2fft(NormalizedPrint,gabor,'same');
                        [disk,vector]=sector_norm(ComponentPrint,1,0);
                        finger_code{angle+1}=vector(1:n_sectors);
                        vettore_in(angle*n_sectors+1:(angle+1)*n_sectors)=finger_code{angle+1};
                    end

                    
                    vettore_a=zeros(num_disk*n_sectors,1);
                    vettore_b=zeros(num_disk*n_sectors,1);
                    best_matching=zeros(fp_number,1);
                    valori_rotazione=zeros(n_arcs,1);
                   
                    for scanning=1:fp_number
                        fcode1=data{scanning,1};
                        fcode2=data{scanning,2};
                        for rotazione=0:(n_arcs-1)
                            p1=fcode1;
                            p2=fcode2;
                            
                            for conta_disco=1:num_disk
                                disco1=p1{conta_disco};
                                disco2=p2{conta_disco};
                                for old_pos=1:n_arcs
                                    new_pos=mod(old_pos+rotazione,n_arcs);
                                    if new_pos==0
                                        new_pos=n_arcs;
                                    end
                                    for conta_bande=0:1:(n_bands-1)
                                        disco1r(new_pos+conta_bande*n_arcs)=disco1(old_pos+conta_bande*n_arcs);
                                        disco2r(new_pos+conta_bande*n_arcs)=disco2(old_pos+conta_bande*n_arcs);
                                    end
                                end
                                p1{conta_disco}=disco1r;
                                p2{conta_disco}=disco2r;
                            end
                            
                            for old_disk=1:num_disk
                                new_disk=mod(old_disk+rotazione,num_disk);
                                if new_disk==0
                                    new_disk=num_disk;
                                end
                                pos=old_disk-1;
                                vettore_a(pos*n_sectors+1:(pos+1)*n_sectors)=p1{new_disk};
                                vettore_b(pos*n_sectors+1:(pos+1)*n_sectors)=p2{new_disk};
                            end
                            d1=norm(vettore_a-vettore_in);
                            d2=norm(vettore_b-vettore_in);
                            if d1<d2
                                val_minimo=d1;
                            else
                                val_minimo=d2;
                            end
                            valori_rotazione(rotazione+1)=val_minimo;
                        end
                        [minimo,posizione_minimo]=min(valori_rotazione);
                        best_matching(scanning)=minimo;
                    end
                    [distanza_minima,posizione_minimo]=min(best_matching);
                    results_img(scan_core) = posizione_minimo;
                    results_dis(scan_core) = distanza_minima;
                    
                end
                
                [distanza_minima,posizione_minimo] = min(results_dis);
                dito_minimo                        = results_img(posizione_minimo);

                message=strcat('Recognized fingerprint:',namefile_vector{dito_minimo});
               
                message=strcat('Location:',path_vector{dito_minimo});
               


                message=strcat('The Nearest Fingerprint Is  : ',num2str(dito_minimo),...
                    ' With a distance Of : ',num2str(distanza_minima));
               
                 name=namefile_vector{dito_minimo};
                 pathx=path_vector{dito_minimo};
                 resultx=strcat(pathx,name);
                 imgx=imread(resultx);
                 
                  
            imshow(imgx);
                msgbox(message,'DataBase Info','help');
            end
        else
            message='DataBase is empty. No check is possible.';
            msgbox(message,'FingerCode DataBase Error','warn');
        end



function pushbutton4_Callback(hObject, eventdata, handles)

if (exist('FingerPrint_database.dat')==2)
            button = questdlg('Do you want to Delete Database?');
            if strcmp(button,'Yes')
                delete('FingerPrint_database.dat');
                msgbox('Database was succesfully Removed.','Database removed','help');
            end
        else
            warndlg('Database is empty.',' Warning ')
        end

 
set(handles.edit3,'String',0);




function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
