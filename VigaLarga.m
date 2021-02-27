%Programa para simular una viga larga

LongitudViga    = input('Introduzca la longitud de la viga :');
AltitudViga     = input('Introduzca la altitud de la viga :');
Bloques         = input('Introduzca la cantidad de celdas :');

%Deje predefinido un modulo de Young y un area transverzal

AreaTransversal = 0.01;
E               = 1.5e11;
kVertical       = (AreaTransversal*E)/AltitudViga;
kHorizontal     = (AreaTransversal*E)/(LongitudViga/Bloques);
kDiagonal       = (AreaTransversal*E)/(sqrt((AltitudViga)^2+(LongitudViga/Bloques)^2));
anguloDiagonal  = atan(AltitudViga/(LongitudViga/Bloques));

%Creo la matriz Global y la de Elementos 

MatrizGlobal = MatrizGlobal(Bloques*2+2,2);

%Creo todos los Elementos verticales:

for i=1:Bloques+1
    if i==1
        El = ELEMENTO(1,3,kVertical,pi/2);
    elseif i==Bloques+1
        El = [El ELEMENTO(2,Bloques+3,kVertical,pi/2)];
    else
        El = [El ELEMENTO(i+2,i+2+Bloques,kVertical,pi/2)];
    end
end

%Creo todos los Elementos horizontales:

for i=1:Bloques*2
    if i<=Bloques
        El = [El ELEMENTO(i+2,i+3,kHorizontal,0)]; %Elementos de arriba
    elseif (i==Bloques+1)
        El = [El ELEMENTO(1,Bloques+4,kHorizontal,0)]; %Elemento de abajo a la izquierda
    elseif (i==Bloques*2)
        El = [El ELEMENTO(2,Bloques*2+2,kHorizontal,0)]; %Elemento de abajo a la derecha
    else
        El = [El ELEMENTO(i+2,i+3,kHorizontal,0)]; %Elementos de abajo no esquinados
    end
end

%Creo todos los Elementos a "45°":

for i=1:Bloques
    if i==1
        El = [El ELEMENTO(1,4,kDiagonal,anguloDiagonal)]; %Elementos de arriba
    else
        El = [El ELEMENTO(Bloques+2+i,i+3,kDiagonal,anguloDiagonal)]; %Elementos de abajo no esquinados
    end
end

%Creo todos los Elementos a "-45°":

for i=1:Bloques
    if i==Bloques
        El = [El ELEMENTO(Bloques+2,2,kDiagonal,-anguloDiagonal)]; %Elementos de arriba
    else
        El = [El ELEMENTO(Bloques+3+i,i+2,kDiagonal,-anguloDiagonal)]; %Elementos de abajo no esquinados
    end
end

%Creo Los nodos para cada punto de la viga que se discretizó

varx  = [];vary = [];%para el grafico de nodos
Nodos = [];
for i=1:Bloques*2+2
    if(i==1)
        Nodos = [Nodos NODO2D(0,0)];
    
    elseif (i==2)
        Nodos = [Nodos NODO2D(LongitudViga,0)];
    
    elseif (i>=3 && i<=3+Bloques)
        Nodos = [Nodos NODO2D((i-3)*(LongitudViga/Bloques),AltitudViga)];
        
    else
        Nodos = [Nodos NODO2D((i-(3+Bloques))*(LongitudViga/Bloques),0)];
    end
    varx = [varx Nodos(i).x]; %para el grafico de nodos
    vary = [vary Nodos(i).y]; %para el grafico de nodos
end



%Ensamblo la matriz global 

[s1,s2] = size(El);
for i=1:s2
    MatrizGlobal=insercion(MatrizGlobal,El(i).matrizLocal,El(i).nodos(1),El(i).nodos(2));
end

%Armo la submatriz KF

s3 = size(MatrizGlobal);
KF = MatrizGlobal(5:s3(1),5:s3(1));

%Armo la submatriz KEF

KEF = MatrizGlobal(1:4,5:s3(1));

%Creo la matriz columna de fuerzas conocidas

F      = [];
fuerza = input('Por favor introduzca la fuerza soportada :');

% -------------Grafico para ubicar visualmente los nodos ------------
figure(1)
plot(varx,vary,"color","none","marker","o","markerFaceColor","k"); %Grafico para identificar el nodo
hold on
for i=1:Bloques*2+2
    text(Nodos(i).x,Nodos(i).y,['\leftarrow' num2str(i)]);
end
xlim([-LongitudViga/25 LongitudViga+LongitudViga/25])
ylim([-AltitudViga/5 AltitudViga+AltitudViga/5])
hold off
%------------------------------------------------------------------

choice = menu('Elija una forma de carga','Puntual','Distribuida');

if (choice == 1)
    Nodo   = input('Por favor introduzca el nodo que la soporta :');
    dim    = input('Por favor introduzca la dimension en la que esta la fuerza siendo 1:x y 2:y :');

    for i=1:((Bloques*2+2)*2)-4
        if (i==Nodo*2-6+dim) 
            F = [F;fuerza];
        else
            F = [F;0];
        end
    end
    
else
    carga = fuerza/(Bloques+1);
    
    for i=1:((Bloques*2+2)*2)-4   
        if (mod(i,2) == 0 && i <= (Bloques+1)*2) 
            F = [F;carga];
        else
            F = [F;0];
        end
    end
end

%Creo el sistema y lo resuelvo para posiciones desconocidas
sistema    = [KF F];
sistema    = rref(sistema);
tam        = size(sistema);
posiciones = sistema(:,tam(2));

%Luego utilizo las posiciones encontradas guardadas en posiciones para usarla
%para encontrar las reacciones con la submatriz KEF
reacciones = KEF * posiciones;

for i=2:2:4
    fprintf("Las reacciones para el nodo %d es:(%d,%d) \n",i/2,reacciones(i-1),reacciones(i));
end

for i=2:2:(Bloques*2+2)*2-4
    fprintf("Las posiciones respecto del punto de equilibro para el nodo %d son:(%d,%d) \n",i/2+2,posiciones(i-1),posiciones(i));
    Nodos((i/2)+2).x   = Nodos((i/2)+2).x + posiciones(i-1);
    Nodos((i/2)+2).eqx =  posiciones(i-1);
    Nodos((i/2)+2).y = Nodos((i/2)+2).y + posiciones(i);
    Nodos((i/2)+2).eqy =  posiciones(i);
    
end

%Calculo los esfuerzos para los elementos


for i=1:s2
    %Uso la funcion esfuerzos, poniendo en ella los nodos relacionados
    %y el elemento.
    El(i).esfuerzo = esfuerzos(Nodos(El(i).nodos(1)) , Nodos(El(i).nodos(2)) , El(i));
    
end



%Seccion de ploteo del sistema (En construccion) 
matCx = [0;0]; matCy = [0;0];
matTx = [0;0]; matTy = [0;0];
for i=1:s2
    if(El(i).esfuerzo < 0)
        a     = [Nodos(El(i).nodos(1)).x ; Nodos(El(i).nodos(2)).x];
        matCx = [matCx a];
        a     = [Nodos(El(i).nodos(1)).y ; Nodos(El(i).nodos(2)).y];
        matCy = [matCy a];
    else
        a     = [Nodos(El(i).nodos(1)).x ; Nodos(El(i).nodos(2)).x];
        matTx = [matTx , a];
        a     = [Nodos(El(i).nodos(1)).y ; Nodos(El(i).nodos(2)).y];
        matTy = [matTy , a];
    end
end

figure(2)
plot(matCx,matCy,'r',matTx,matTy,'b');
hold on
plot(matCx,matCy,'color','none','marker','o','markerFaceColor','k','markerSize',5);
plot(matTx,matTy,'color','none','marker','o','markerFaceColor','k','markerSize',5);


for i=1:Bloques*2+2
    text(Nodos(i).x,Nodos(i).y,['\leftarrow' num2str(i)]);
end
xlim([-LongitudViga/25 LongitudViga+LongitudViga/25])
ylim([-AltitudViga/2 AltitudViga+AltitudViga/2])

hold off