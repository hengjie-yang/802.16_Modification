function fraction = Compute_transfer_function(NetFile, Start, Stop)


[Num, Den] = mason(NetFile, Start, Stop);
Num = str2sym(Num);
Den = str2sym(Den);

disp(['Start node: ',num2str(Start),' End node: ',num2str(Stop),' Transfer function: ']);
fraction = simplifyFraction(Num/Den);