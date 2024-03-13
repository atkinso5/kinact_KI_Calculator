#Code for calculating kinact and KI from an IC50 curve
#Please refer to the following paper: Mader & Keillor, ACS Med. Chem. Lett. 2024, 15, xxxx.
#The Python code was developed and is maintained by Bethany Atkinson at Dunad Therapeutics.
#The code is free to use but please cite "Mader & Keillor, ACS Med. Chem. Lett. **2024**, *15*, xxxx" and the GitHub link when using it.


import numpy as np 
from scipy.optimize import least_squares 

class KineticParameterCalculator:
    def __init__(self, array, AddSub, EnzConc, kcat, Km, PreIncVol, IncVol, kinact=1, KI=1): #Initial estimates of kinact and KI default to 1 but can be overwritten to check for convergence from different starting points 
        self.array = array

        self.AddSub = AddSub
        self.EnzConc = EnzConc
        self.kcat = kcat
        self.Km = Km
        self.PreIncVol = PreIncVol
        self.IncVol = IncVol
        self.DilFact = PreIncVol / IncVol

        self._kinact = kinact # 1
        self._KI = KI # 1

        #call methods to get predicted kinact and KI values 
        self.get_response_coeff()
        self.calculate_kinact_KI()
        self.calculate_Pred_signal([self._kinact, self._KI])
        self.update_array_w_pred_signal()
        self.stats_of_optimisation()

        #print calculated kinact and KI values:
        print("kinact = " + str(self._kinact) + "\n" + "KI = " + str(self._KI))


    # simulates reaction and estimate final product concentration 
    def PreIncEndPoint(self, inputs, PreIncTime, IncTime, InhConc): # inputs = [kinact, KI ] 
        kinact = inputs[0]
        KI = inputs[1]  

        EnzConc = self.EnzConc #Do not want to override as will change value output from PreIncEndPoint when recall the function 

        #Set the granularity of each phase of the simulation to 100  finely-divided time intervals
        dPreTime = PreIncTime / 100
        dIncTime = IncTime / 100

        SubConc = 0 #set initial concentration 
        ProdConc = 0 #set initial concentration

        for i in range(100):    #For each time interval dPreTime (mimics a simulation over the PreIncTime with EnzConc and InhConc increasing each 100th of the simulation)
            #First calculate instantaneous rate of E-I formation
            EIRate = kinact * InhConc / (InhConc + KI) * EnzConc
                
            #Now calculate incremental changes in concentrations, multiplying rates by time interval (dPreTime)
            dEIConc = EIRate * dPreTime     #E and I both decrease by the same amount
            if dEIConc > EnzConc:
                dEIConc = EnzConc         #This protects from EnzConc going below zero

            if dEIConc > InhConc:
                dEIConc = InhConc         #This protects from InhConc going below zero

            #Then calculate new concentrations, at the end of this time interval:
            EnzConc = EnzConc - dEIConc  #Conc decrease for enzyme.
            InhConc = InhConc - dEIConc  #Conc decrease for inhibitor.


        #Now account for addition of substrate and dilution of Enz and Inh
        SubConc = SubConc + self.AddSub 
        EnzConc = EnzConc * self.DilFact
        InhConc = InhConc * self.DilFact


        #Now simulate Incubation phase
        for j in range(100):   #For each time interval dIncTime (mimics a simulation over the IncTime with EnzConc and InhConc increasing each 100th of the simulation)
            
            #First calculate instantaneous rate of prod formation, accounting for competitive inhibition:
            InstRate = self.kcat * EnzConc * SubConc / (SubConc + self.Km * (1 + InhConc / KI))
            #Then calculate instantaneous rate of E-I formation, accounting for competition with substrate:
            EIRate = kinact * (InhConc / (InhConc + KI * (1 + SubConc / self.Km))) * EnzConc
            
            #Now calculate incremental changes in concentrations, multiplying rates by time interval (dIncTime)
            dSPConc = InstRate * dIncTime  #Sub and Prod change by the same amount
            if dSPConc > SubConc:
                dSPConc = SubConc  #This protects from SubConc going below zero
            
            dEIConc = EIRate * dIncTime     #E and I both decrease by the same amount
            if dEIConc > EnzConc:
                dEIConc = EnzConc         #This protects from EnzConc going below zero

            if dEIConc > InhConc:
                dEIConc = InhConc         #This protects from InhConc going below zero


            #Then calculate new concentrations, at the end of this time interval:
            SubConc = SubConc - dSPConc     #Conc decrease for substrate
            ProdConc = ProdConc + dSPConc   #Conc increase for product 
            EnzConc = EnzConc - dEIConc     #Conc decrease for enzyme
            InhConc = InhConc - dEIConc     #Conc decrease for inhibitor

        PreIncEndPoint = ProdConc     #Return final product concentration
        return PreIncEndPoint 


    #estimates the response coefficient to scale pred_signal
    def get_response_coeff(self):
        #PredSignal * coeff = 100 (when [Inhibitor] = 0) 
        #coeff = 100 / PredSignal 
        response_coeff = 100 / self.PreIncEndPoint([1,1], self.array[0,0], self.array[0,1], 0)
        self.response_coeff = response_coeff
        return response_coeff 


    #calculates predicted signal for entire array (applies PreIncEndPoint to each line). Requires input of kinact and KI 
    def calculate_Pred_signal(self, inputs):
        Pred_signal = []
        for i in range(len(self.array)): #better way than looping to apply to full array? May not be possible as can't vectorise the PreIncEndPoint function 
                Pred_signal.append(self.PreIncEndPoint(inputs, self.array[i,0], self.array[i,1], self.array[i,2]) * self.response_coeff)
        self.pred_signal = Pred_signal
        return Pred_signal 


    # return calculated kinact and KI value 
    def calculate_kinact_KI(self):
        inputs = [self._kinact, self._KI]
        def calculate_residuals(inputs): #need function as input to scipy least_squares which takes [kinact, KI] as input and calculates the residuals 
            Pred_signal = self.calculate_Pred_signal(inputs)
            Residuals = self.array[:,3] - Pred_signal 
            return Residuals

        self.least_squares_results = least_squares(calculate_residuals, inputs, bounds = [0, np.inf]) #bounds ensure kinact and KI are >= 0 
        self._kinact = self.least_squares_results.x[0] #override kinact with calculated value 
        self._KI = self.least_squares_results.x[1] #override KI with calculated value
        return self._kinact, self._KI


    #to calculate predicted signal once calculated kinact and KI
    def update_array_w_pred_signal(self):
        Pred_signal = self.calculate_Pred_signal([self._kinact, self._KI])
        self.array = np.concatenate([self.array, np.array(Pred_signal).reshape(len(Pred_signal), 1)], axis=1) #add column of predicted values to numpy array 

    def stats_of_optimisation(self):
        Pred_signal = self.calculate_Pred_signal([self._kinact, self._KI])
        residuals = self.array[:,3] - Pred_signal 
        sum_sq_residuals = np.sum(residuals**2)
        sum_sq_total = np.sum((self.array[:,3] - self.array[:,3].mean())**2) #signal - mean(signal)
        self.R2 = 1 - (sum_sq_residuals / sum_sq_total)
        self.RMS_error = np.sqrt(sum_sq_residuals / len(self.array[:,3]))
