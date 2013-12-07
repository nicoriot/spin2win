function idat = Data_2_idat(Data)
    for i = 1:Data.nShells
        idat(1,i)  = Data.Shells(i).d.i;
        idat(2,i)  = Data.Shells(i).d.o;
        idat(3,i)  = Data.Shells(i).E.c;
        idat(4,i)  = Data.Shells(i).E.r;
        idat(5,i)  = Data.Shells(i).v;
        idat(6,i)  = Data.Shells(i).rho;
        idat(7,i)  = Data.Shells(i).mu;
        idat(9,i)  = Data.Shells(i).G.rz;
        idat(10,i) = Data.Shells(i).Cost;
    end
    idat(8,1)     = Data.nRadial;
    idat(8,2)     = Data.height;
    idat(8,3)     = Data.rot_speed.max;
    idat(8,4)     = Data.rot_speed.min;
    idat(8,5)     = Data.nShells;
end