#############################################################################
#Question 1
 
df=read.table("copula-hw1.txt", header=TRUE)
epsilon=1

#sub sample X, Y into the derivative functions, using the given mu1, mu2 and theta
#return the calcualted value of mu1, mu2 and theta
derivative_density <- function(X, Y, mu1, mu2, theta)
{
	densityX=dnorm(X - mu1)
	densityY=dnorm(Y - mu2)
	probX=pnorm(X - mu1)
	probY=pnorm(Y - mu2)
    
	dF_dmu1_numerator=(2 * theta * exp(-1 * theta * probX) *densityX * (1 - exp(-1 * theta * probY)))    
	dF_dmu1_denominator=(1 - exp(-1 * theta) - (1-exp(-1 * theta * probX)) * (1 - exp(-1 * theta * probY)))
	dF_dmu1=(d_mu1_part1=(X - mu1) + theta * densityX) - dF_dmu1_numerator/dF_dmu1_denominator
    
	dF_dmu2_numerator=(2 * theta * exp(-1 * theta * probY) * densityY * (1 - exp(-1 * theta * probX)))    
	dF_dmu2_denominator=(1 - exp(-1 * theta) - (1 - exp(-1 * theta * probY)) * (1 - exp(-1 * theta * probX)))
	dF_dmu2=((Y - mu2)+theta * densityY) - dF_dmu2_numerator/dF_dmu2_denominator
    
	dF_dtheta_numerator=2 * (exp(-1 * theta)+((-1 * probX * exp(-1 * theta * probX) * (1 - exp(-1 * theta*probY))) + (-1 * probY * exp(-1 * theta * probY) * (1 - exp(-1 * theta * probX)))))
	dF_dtheta_denominator=(1 - exp(-1 * theta) - (1 - exp(-1 * theta * probY)) * (1 - exp(-1 * theta * probX)))    
	dF_dtheta=(1/theta + exp(-1 * theta)/(1 - exp(-1 * theta)) - (probX + probY)) - dF_dtheta_numerator/dF_dtheta_denominator
	return(c(dF_dmu1,dF_dmu2,dF_dtheta))    
}

#the log the density function
log_density <- function(X, Y, mu1, mu2, theta)
{
	densityX=dnorm(X - mu1)
	densityY=dnorm(Y - mu2)
	probX=pnorm(X - mu1)
	probY=pnorm(Y - mu2)
	density_numerator= densityX * densityY * theta * (1 - exp(-1 * theta)) * exp(-1 * theta * (probX + probY))
	density_denominator =(1 - exp(-1 * theta) - (1 - exp(-1 * theta * probX)) * (1 - exp(-1 * theta * probY))) **2
	density=log(density_numerator/density_denominator)
	return(density)
}


derivative_log_likelihood <- function(df, mu1, mu2, theta)
{
    result=c()
    for(i in 1:200)
    {

        temp=derivative_density(df[i,1], df[i,2], mu1, mu2, theta)
        result=rbind(result,(derivative_density(df[i,1], df[i,2], mu1, mu2, theta)))
    }
    
    return(apply(result,2,sum))
}

#second derivative was too hard to calculate, so here I am using the secant method to estimate the Hassian matrix 
#for anything that is not the diagnal of the hessian_matrix, I am using 0, to help estimate due to the complex calculation
hessian_matrix <- function(mu10, mu11, mu20, mu21, theta0, theta1, xi, yi)
{
	dd1=derivative_density(xi, yi, mu10, mu20, theta0)
	dd2=derivative_density(xi, yi, mu11, mu21, theta1)
	
	denominator11=(mu11-mu10)/(dd2[1] - dd1[1])
	denominator22=(mu21-mu20)/(dd2[2] - dd1[2])
	denominator33=(theta1-theta0)/(dd2[3] - dd1[3])
	
	v1=c(1/denominator11,0,0)
	v2=c(0,1/denominator22,0)
	v3=c(0,0,1/denominator33)
	result=rbind(v1,v2,v3)
   	return(result)
}


inverse_hessian <- function(df, mu10, mu11, mu20, mu21, theta0, theta1)
{
    result=matrix(rep(0,9),ncol=3,nrow=3)
    for(i in 1:200)
    {
    	non_inverse_result=hessian_matrix(mu10, mu11, mu20, mu21, theta0, theta1, df[i,1], df[i,2])
        result=result+non_inverse_result
    }
   
    return(solve(result))
}


#select random values to start the optimization
result_temp =c(0.5,0.5,0.5)
result =c(1,2,3)

likelihood =derivative_log_likelihood(df, 1, -1, 1)

#here is the optimization loop, where it will stop once error is smaller than 0.00001

for(k in 1:1000)
{
    mu10= result_temp[1]
    mu11= result[1]
    
    mu20= result_temp[2]
    mu21= result[2]    
    
    theta0= result_temp[3]
    theta1= result[3]
    
    hessian=inverse_hessian(df, mu10, mu11, mu20, mu21, theta0, theta1)

    result_new= result - 0.1* hessian %*% likelihood
    mu1= result_new[1]
    mu2= result_new[2]
    theta= result_new[3]
    
    result_temp = result
    likelihood=derivative_log_likelihood(df, mu1, mu2, theta)
    epsilon=max((abs((result[1] - result_new[1])/result[1])),(abs((result[2] - result_new[2])/result[2])),(abs((result[3] - result_new[3])/result[3])))
    result = result_new
    if(epsilon<0.00001){break}
}
message("Q1")
message("mu1 is: ", mu1)
message("mu2 is: ", mu2)
message("theta is: ", theta)


#############################################################################
#Question 2
antithetic_boot = read.table("antithetic-boot.txt")
#randomize the original vector
original <- as.numeric(readLines("antithetic-boot.txt"))
rand1 <- as.numeric(antithetic_boot[sample(nrow(antithetic_boot), replace=FALSE),])
rand2 <- as.numeric(antithetic_boot[sample(nrow(antithetic_boot),replace=FALSE),])

#calculate and create new columns with the first 3 columns
df = data.frame(original,rand1,rand2)
df$YjYr1j <- df$original * df$rand1
df$YjYr2j <- df$original * df$rand2
df$Yr1jYr2j <- df$rand1 * df$rand2

#find sum function
sum_all <- function(df)
{
	df$row_sum <- df$YjYr1j + df$YjYr2j + df$Yr1jYr2j
	total_sum_result = sum(df$row_sum)
	return(total_sum_result)
}

#find possible rows to exchange, this is defined as my local search
exchange_rand1 <- function (df, row1, row2)
{
	temp_df = df
	temp_df$rand1[row1] = df$rand1[row2]
	temp_df$rand1[row2] = df$rand1[row1]
	temp_df$YjYr1j <- temp_df $original * temp_df$rand1
	temp_df$YjYr2j <- temp_df $original * temp_df$rand2
	temp_df$Yr1jYr2j <- temp_df $rand1 * temp_df$rand2
	return(temp_df)
}
exchange_rand2 <- function (df, row1, row2)
{
	temp_df = df
	temp_df$rand2[row1] = df$rand2[row2]
	temp_df$rand2[row2] = df$rand2[row1]
	temp_df$YjYr1j <- temp_df $original * temp_df$rand1
	temp_df$YjYr2j <- temp_df $original * temp_df$rand2
	temp_df$Yr1jYr2j <- temp_df $rand1 * temp_df$rand2
	return(temp_df)
}
opt <- function(df)
{
	#initialize sum
	current_sum  = sum_all(df)
	new_df = df
	iteration = 0
	error = 1
	len = 64
	
	while(iteration < 100)
	{
		#sum before changes
		init_sum = sum_all(new_df)
		#for column rand1, if find a swap that will optimize the result, switch
		for(j in 1:(len-1))
		{
			for(k in (j+1):len)
			{

				temp_df = exchange_rand1(new_df, j, k)
				temp_sum = sum_all(temp_df) 
				if(temp_sum < current_sum){
					new_df = temp_df
					current_sum = temp_sum
				}
			}
		}
		
		#for column rand2, if find a swap that will optimize the result, switch		
		for(j in 1:(len-1))
		{
			for(k in (j+1):len)
			{
				temp_df = exchange_rand2(new_df, j, k)
				temp_sum = sum_all(temp_df) 
				if(temp_sum < current_sum){
					new_df = temp_df
					current_sum = temp_sum
				}				
			}
		}
		iteration = iteration + 1
		error = abs(init_sum - current_sum)/abs(current_sum)
		#print(error)
		#print(iteration)
		#stop interation once the error is small enough
		if(error<0.00001){break}
	}
	return(new_df)
}

new_df = opt(df)
sum_all(new_df)
message("Q2, sum is: ", sum_all(new_df))
message("To check the optimized columns (rand1 and rand2), type new_df and hit enter")


##########################################################################



#Question 3

distance <- function(row, center)
{
	diff = sqrt(sum((row-center)**2))
	return(diff)	
}
sum_sqr <- function(df)
{
	data = list(df[,2:14])
	len = dim(df)[1]
	center= aggregate(data, by=list(df$h),FUN='mean', na.rm=TRUE)
	temp_df = df
	for(i in 1:len)
	{
		group = temp_df[i,"h"]
		group_center = center[center$Group.1==group,2:14]
		distance = distance(temp_df[i,2:14],group_center)
		temp_df$sum_sqr[i] = distance**2
	}
	result=sum(temp_df$sum_sqr)
	return (result)
	
}


wine = read.table("wine.dat", header = TRUE)
#create a new column giving out random num from 1-3 to create random groups of 3
h <- as.numeric(sample(1:3, 178,replace=TRUE))
df = data.frame(wine)
df$h = h
len = 178
h1 <- subset(df, h == 1, 2:15)
h2 <- subset(df, h == 2, 2:15)
h3 <- subset(df, h == 3, 2:15)

data = list(df[,2:14])
current_center= aggregate(data, by=list(df$h),FUN='mean', na.rm=TRUE)

new_center = current_center

for(i in 1:100)
{
	h1_center=new_center[new_center $Group.1==1,2:14]
	h2_center=new_center[new_center $Group.1==2,2:14]
	h3_center=new_center[new_center $Group.1==3,2:14]
	
	current_sum_sqr = sum_sqr(df)
	
	#find distance of each point to the 3 group centers, 
	#then put it in the group that is the closest
	#this way we would be able to diverge into the best possible groups
	for(j in 1:len)
	{
		j_data = df[j, 2:14]
		d1 = distance(j_data, h1_center)
		d2 = distance(j_data, h2_center)
		d3 = distance(j_data, h3_center)
		df$d1[j] = d1
		df$d2[j] = d2
		df$d3[j] = d3
		
		df$h[j] = which.min(c(d1,d2,d3))
	}
	
	new_center= aggregate(data, by=list(df$h),FUN='mean', na.rm=TRUE)
	new_sum_sqr = sum_sqr(df)
	#message("i: ", i)
	error = abs((new_sum_sqr - current_sum_sqr) / current_sum_sqr)
	#print(error)
	if (error < 0.0001){break}
	
}
message("Q3, sum of square is: ", new_sum_sqr)
message("To check the new groups, type df and hit enter, column h is the optimized grouping calculated")

