#-------------------Valores de Entrada-----------------------

# Operacion con pivoteo parcial
A=rbind(c(0.003000,59.14),c(5.291,-4))
b=c(59.17,46.78)

#A=rbind(c(3,1,1),c(1,-4,1),c(1,1,-3))
#b=c(1,2,1)

#A = rbind(c(1,1,0,3), c(2,1,-1,1),c(3,-1,-1,2), c(-1,2,3,-1))
#b = c(4,1,-3,4)

#-------------------Matriz fuerte en sentido diagoanl-----------------------

Dominant<-function(A)
{
	# Teorema. 11. Una matriz A estrictamente dominante en sentido diagonal 
	# es no singular. Más aún, en est caso podemos realizar la eliminación 
	# gaussiana de cualquier sistema lineal de la forma Ax = b para obtener 
	# su solución única sin intercambios de filas ni columnas, y los cálculos 
	# son estables respecto al crecimiento de los errores de redondeo.
	
	print('La Matriz es dominante?')
	fila = dim(A)[1]							# Numero de fila de la matriz A
	columna = dim(A)[2]							# Numero de columna de la matriz A
	if(fila == columna)							# Es Cuadrada
	{	
		i=1
		while(i<=fila && (abs(A[i,i])>sum(abs(A[i,i:fila]))-abs(A[i,i]))) i=i+1
		if(i>fila) return("La matriz es dominante en sentido diagonal")
		if(i<=fila) return("La matriz no es dominante en sentido diagonal")
	}
	return('La matriz no es cuadrada')
}

#------------Eliminacion gaussiana hacia adelante con sustitucion hacia atras--------------------

elimGauss_Adelante_SustAtras<-function(A,b)
{
	print('Eliminacion gaussiana hacia adelante con sustitucion hacia atras')
	if (det(A)==0) return('Sistema singular')   # caso imposible de n−1 ecuaciones con n incógnitas.
	fila = dim(A)[1]							# Numero de fila de la matriz A
	columna = dim(A)[2]							# Numero de columna de la matriz A
	if(fila == columna)							# Es Cuadrada
	{	
		MA=cbind(A,b)							# Matriz Ampliada
		print(MA)								# Imprimir la Matriz Ampliada
		numeroOperaciones(A)					# Imprime el numero de operaciones a realizar
		na=dim(MA)[2]							# Numero de columna de la matriz Ampliada
		for (k in 1:(fila-1))					# Repetir Miestra 
			for (i in (k+1):fila)				
 			{ 
 				if(MA[k,k]==0) 
 					return('Division entre cero')# Para visualizar el metodo falla
				factor=MA[i,k]/MA[k,k]
 				MA[i,k:na]=MA[i,k:na]-factor*MA[k,k:na]
		  	}
		x=numeric(fila)
		x[fila]=MA[fila,na]/MA[fila,fila]
		for ( i in (fila-1):1)
			x[i]=(MA[i,na]-sum(MA[i,(i+1):fila]*x[(i+1):fila]))/MA[i,i] 
		return(x)
	}
	return('La matriz no es cuadrada')
}

#------------Eliminacion Gaussiana hacia atras sustitucion hacia adelante--------------------

elimGauss_Atras_SustAdelante <- function(A,b)
{
	print('Eliminacion Gaussiana hacia atras sustitucion hacia adelante')
	if (det(A)==0) return('Sistema singular')   # caso imposible de n−1 ecuaciones con n incógnitas.
	fila = dim(A)[1]							# Numero de fila de la matriz A
	columna = dim(A)[2]							# Numero de columna de la matriz A
	if(fila == columna)							# Es Cuadrada
	{
		MA=cbind(A,b)							# Matriz Ampliada
		print(MA)								# Imprimir la Matriz Ampliada
		numeroOperaciones(A)					# Imprime el numero de operaciones a realizar
		na=dim(MA)[2]							# Numero de columna de la matriz Ampliada
		for(k in fila:2)
		{
			for(i in (k-1):1)
			{
				if(MA[k,k]==0) 
 					return('Division entre cero')# Para visualizar el metodo falla
				factor=MA[i,k]/MA[k,k]
				MA[i,na:1]=MA[i,na:1]-factor*MA[k,na:1]	
			}
		}
		x=numeric(fila)
		x[1]=MA[1,na]/MA[1,1]
		for(i in 2:fila)
			x[i]=(MA[i,na]-sum(MA[i,1:(i-1)]*x[1:(i-1)]))/MA[i,i]
		return(x)
	}	
	return('La matriz no es cuadrada')
}

#------------Pivoteo parcial hacia atras------------------

pivoteoParcialAtras<- function(MA,n,k)
{
	mayor=which.max(abs(MA[1:k,k]))
	if(MA[mayor,k]==0)
	{
		print("El sistema no tiene solucion unica")
		return()
	}
	if((mayor) != k)
	{
		print(c('Fila a intercambiar: ',MA[mayor,]))
		pivo=MA[k,]
		MA[k,]=MA[mayor,]						# Intercambiamos Fila k por mayor
		MA[mayor,]=pivo 						# Intercambiamos Fila
	}
	print('Matriz resultante')
	print(MA)
	return (MA)
}

#-------- ElimGaus hacia atras con pivoteo Parcial----------

elimGPivoteoParcialAtras <- function(A,b)
{
	print('ElimGaus hacia atras con pivoteo Parcial')
	if (det(A)==0) return('Sistema singular')   # caso imposible de n−1 ecuaciones con n incógnitas.
	fila = dim(A)[1]							# Numero de fila de la matriz A
	columna = dim(A)[2]							# Numero de columna de la matriz A
	if(fila == columna)							# Es Cuadrada
	{
		MA=cbind(A,b)							# Matriz Ampliada
		print(MA)								# Imprimir la Matriz Ampliada
		na=dim(MA)[2]							# Numero de columna de la matriz Ampliada
		for(k in fila:2)
		{
	 		MA <- pivoteoParcialAtras(MA,fila,k)
			for(i in (k-1):1)
			{
				if(MA[k,k]==0) 
 					return('Division entre cero')# Para visualizar el metodo falla
				factor = MA[i,k]/MA[k,k]
				MA[i,na:1]=MA[i,na:1]-factor*MA[k,na:1]	
			}
		}
		x=numeric(fila)
		x[1]=MA[1,na]/MA[1,1]
		for(i in 2:fila)
			x[i]=(MA[i,na]-sum(MA[i,1:(i-1)]*x[1:(i-1)]))/MA[i,i]
		return(x)
	}
	return('La matriz no es cuadrada')	
}	

#------------Pivoteo Parcial Hacia adelante------------------

pivoteoParcialAdelante<- function(MA,n,k)
{
	mayor=which.max(abs(MA[k:n,k]))
	if(MA[mayor+k-1,k]==0)
	{
		print("El sistema no tiene solucion unica")
		return()
	}
	if((mayor+k-1) != k)
	{
		print(MA[mayor+k-1,])
		pivote=MA[k,]
		MA[k,]=MA[mayor+k-1,]
		MA[mayor+k-1,]=pivote
	}
	return (MA)
}

#-------- ElimGaus Hacia adelante con Pivoteo Parcial----------

elimGPivoteoParcialAdelante <- function(A,b)
{
	print('ElimGaus Hacia adelante con Pivoteo Parcial')
	if (det(A)==0) return('Sistema singular')   # caso imposible de n−1 ecuaciones con n incógnitas.
	fila = dim(A)[1]							# Numero de fila de la matriz A
	columna = dim(A)[2]							# Numero de columna de la matriz A
	if(fila == columna)							# Es Cuadrada
	{
		MA=cbind(A,b)							# Matriz Ampliada
		print(MA)								# Imprimir la Matriz Ampliada
		na=dim(MA)[2]							# Numero de columna de la matriz Ampliada
		for(k in 1:(fila-1))
		{
			MA=pivoteoParcialAdelante(MA,fila,k)
			for(i in (k+1):fila)
			{
				if(MA[k,k]==0) 
 					return('Division entre cero')# Para visualizar el metodo falla
				factor = MA[i,k]/MA[k,k]
				MA[i,k:na]=MA[i,k:na]-factor*MA[k,k:na]		
			}
		}
		x=numeric(fila)
		x[fila]=MA[fila,na]/MA[fila,fila]
		for(i in (fila-1):1)
			x[i]=(MA[i,na]-sum(MA[i,(i+1):fila]*x[(i+1):fila]))/MA[i,i]
		return (x)
	}
	return('La matriz no es cuadrada')	
}

#--------------- Pivoteo escalonado mayor--------------

pivoteoEscaladoMayor <- function(MA,n)
{
	s=numeric(n)
	for(i in 1:n)
	{
		s[i]=max(abs(MA[i,i:n]))
		if(s[i] == 0)
		{
			print("El sistema no tiene solucion unica")
			return ()
		}
	}
	return (s)
}

#--------------- Pivoteo Escalonado hacia adelante---------------

pivoteoEscalonadoAdelante <- function(MA,n,k,s)
{
	mayor=which.max(abs(MA[k:n,k]/s[k:n]))
	if((mayor+k-1) != k)
	{
		esca=MA[k,]
		MA[k,]=MA[mayor+k-1,]
		MA[mayor+k-1,]=esca
	}
	return (MA)
}

#-------- ElimGaus Hacia Adelante con Pivoteo Escalonado----------

elimGPivoteoEscalonadolAdelante <- function(A,b)
{
	print('ElimGaus Hacia Adelante con Pivoteo Escalonado')
	if (det(A)==0) return('Sistema singular')   # caso imposible de n−1 ecuaciones con n incógnitas.
	fila = dim(A)[1]							# Numero de fila de la matriz A
	columna = dim(A)[2]							# Numero de columna de la matriz A
	if(fila == columna)							# Es Cuadrada
	{	
		MA=cbind(A,b)							# Matriz Ampliada
		print(MA)								# Imprimir la Matriz Ampliada
		na=dim(MA)[2]							# Numero de columna de la matriz Ampliada
		s=pivoteoEscaladoMayor(MA,fila)			# Se Obtiene el pivote mayor
		print(s)
		for(k in 1:(fila-1))
		{
		 	MA=pivoteoEscalonadoAdelante(MA,fila,k,s)
			for(i in (k+1):fila)
			{
				if(MA[k,k]==0) 
 					return('Division entre cero')# Para visualizar el metodo falla
				factor=MA[i,k]/MA[k,k]
				MA[i,k:na]=MA[i,k:na]-factor*MA[k,k:na]
			}
		}
		x=numeric(fila)
		x[fila]=MA[fila,na]/MA[fila,fila]
		for(i in (fila-1):1)
			x[i]=(MA[i,na]-sum(MA[i,(i+1):fila]*x[(i+1):fila]))/MA[i,i]
		return (x)
	}
	return('La matriz no es cuadrada')	
}

#--------------- Pivoteo Escalonado hacia atras---------------

pivoteoEscalonadoAtras <- function(MA,n,k,s)
{
	mayor=which.max(abs(MA[1:k,k]/s[1:k]))
	if((mayor) != k)
	{
		esca=MA[k,]
		print('POLICE problema con el indice A22=0')
		print(esca)
		print(mayor)
		print(k)
		MA[k,]=MA[mayor+k-1,] ##MA[k,]=MA[mayor+k-2,] solucion
		MA[mayor+k-1,]=esca
	}
	return (MA)
}

#-------- ElimGaus Hacia atras con Pivoteo Escalonado----------

elimGPivoteoEscalonadolAtras <- function(A,b)
{
	print('ElimGaus Hacia Adelante con Pivoteo Escalonado')
	if (det(A)==0) return('Sistema singular')   # caso imposible de n−1 ecuaciones con n incógnitas.
	fila = dim(A)[1]							# Numero de fila de la matriz A
	columna = dim(A)[2]							# Numero de columna de la matriz A
	if(fila == columna)							# Es Cuadrada
	{	
		MA=cbind(A,b)							# Matriz Ampliada
		print(MA)								# Imprimir la Matriz Ampliada
		na=dim(MA)[2]							# Numero de columna de la matriz Ampliada
		s=pivoteoEscaladoMayor(MA,fila)			# Se Obtiene el pivote mayor
		print(s)
		for(k in fila:2)
		{
		 	MA=pivoteoEscalonadoAtras(MA,fila,k,s)
			for(i in (k-1):1)
			{
				if(MA[k,k]==0) 
 					return('Division entre cero')# Para visualizar el metodo falla
				factor = MA[i,k]/MA[k,k]
				MA[i,na:1]=MA[i,na:1]-factor*MA[k,na:1]	
			}
		}
		x=numeric(fila)
		x[1]=MA[1,na]/MA[1,1]
		for(i in 2:fila)
			x[i]=(MA[i,na]-sum(MA[i,1:(i-1)]*x[1:(i-1)]))/MA[i,i]
		return(x)
	}
	return('La matriz no es cuadrada')	
}

jacobi <- function(A,b,tol)
{
	print(A)
	n=dim(A)[1]
	Q=diag(diag(A),n,n)
	print(Q)
	R=(A-Q)
	x=numeric(n)
	iter=1
	error=tol + 1
	print(round(x,4))
	while(error > tol)
	{
		ant=x
		for(i in 1:n)
		{
			x[i]=(b[i]-sum(R[i, 1:n]*ant[1:n]))/Q[i,i]
			print(round(x,4))
		}
		iter=iter + 1
		error=round(sqrt(sum((x[1:n]-ant[1:n])^2)),3)	               
		print(error)
	}
}	

#------------Conteo de operaciones eliminacion gaussiana--------------------

numeroOperaciones <- function(A){
	k=1
	n=dim(A)[1]
	op=sum(k:n)+sum((k:(n-1))*((k+2):(n+1))+(0:(n-2)))
	print(c('Numero de operaciones:',op))
}

#------------Metodos a ejecutar--------------------

print(Dominant(A))								# Muestra si la matriz es dominante en sentido diagonal

x=elimGauss_Adelante_SustAtras(A,b)
print('Resultado')
print(x)

x=elimGauss_Atras_SustAdelante(A,b)
print('Resultado')
print(x)

x=elimGPivoteoParcialAtras(A,b)
print('Resultado')
print(x)

x=elimGPivoteoParcialAdelante(A,b)
print('Resultado')
print(x)

x=elimGPivoteoEscalonadolAdelante(A,b)
print('Resultado')
print(x)

x=elimGPivoteoEscalonadolAtras(A,b)
print('Resultado')
print(x)

