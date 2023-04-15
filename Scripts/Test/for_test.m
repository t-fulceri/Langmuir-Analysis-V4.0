    Nx = 23
    Ny = 18
    Nr = 5

    N = int16(Nx * Ny * Nr)
    for n = 0:N-1
        nr = int16(mod(n,Nr)+1);
        aux = int16(idivide(n,Nr,'floor'));
        ny = int16(mod(aux,Ny)+1);
        nx = int16(idivide(aux,Ny,'floor')+1);
        
        fprintf(['Total_index = ',num2str(n+1),'\n']);
        fprintf(['X_index = ',num2str(nx),'\n']);
        fprintf(['Y_index = ',num2str(ny),'\n']);
        fprintf(['Repetition_index = ',num2str(nr),'\n\n']);
        
    end