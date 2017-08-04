from libsbml import *
from peitho.errors_and_parsers.abc_sysbio.abcsysbio.relations import *
import os
import re
from peitho.errors_and_parsers.abc_sysbio.abcsysbio_parser.Writer import Writer

class SdeCUDAWriter(Writer):
    def __init__(self, sbmlFileName, modelName="", inputPath="", outputPath=""):
        Writer.__init__(self, sbmlFileName, modelName, inputPath, outputPath)
        self.out_file=open(os.path.join(outputPath,self.parsedModel.name+".cu"),"w")
        
    def mathMLConditionParserCuda(self, mathMLstring):
        """
        Replaces and and or with and_ and or_ in a MathML string.
        Returns the string with and and or replaced by and_ and or_
    
        ***** args *****
    
        mathMLstring:
    
                A mathMLstring
    
        """
        
        andString = re.compile("and")
        orString = re.compile("or")
        mathMLstring = andString.sub("and_", mathMLstring)
        mathMLstring = orString.sub("or_", mathMLstring)

        return mathMLstring
    
    def write(self, method):     
        p=re.compile('\s')
    
        #Write number of parameters and species
        self.out_file.write("#define NSPECIES " + str(self.parsedModel.numSpecies) + "\n")
        self.out_file.write("#define NPARAM " + str(self.parsedModel.numGlobalParameters) + "\n")
        self.out_file.write("#define NREACT " + str(self.parsedModel.numReactions) + "\n")
        self.out_file.write("\n")
    
        #The user-defined functions used in the model must be written in the file
        self.out_file.write("//Code for texture memory\n")
    
        numEvents = len(self.parsedModel.listOfEvents)
        numRules = len(self.parsedModel.listOfRules)
        num = numEvents+numRules
        if num>0:
            self.out_file.write("#define leq(a,b) a<=b\n")
            self.out_file.write("#define neq(a,b) a!=b\n")
            self.out_file.write("#define geq(a,b) a>=b\n")
            self.out_file.write("#define lt(a,b) a<b\n")
            self.out_file.write("#define gt(a,b) a>b\n")
            self.out_file.write("#define eq(a,b) a==b\n")
            self.out_file.write("#define and_(a,b) a&&b\n")
            self.out_file.write("#define or_(a,b) a||b\n")
    
        for i in range(0,len(self.parsedModel.listOfFunctions)):
            self.out_file.write("__device__ float "+self.parsedModel.listOfFunctions[i].getId()+"(")
            for j in range(0, self.parsedModel.listOfFunctions[i].getNumArguments()):
                self.out_file.write("float "+self.parsedModel.functionArgument[i][j])
                if(j<( self.parsedModel.listOfFunctions[i].getNumArguments()-1)):
                    self.out_file.write(",")
            self.out_file.write("){\n    return ")
            self.out_file.write(self.parsedModel.functionBody[i])
            self.out_file.write(";\n}\n")
            self.out_file.write("\n")
    
        self.out_file.write("\n")
    
        self.out_file.write("__device__ void step(float *y, float t, unsigned *rngRegs, int tid){\n")
    
        #write rules and events
        for i in range(0,len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isRate() == True:
                self.out_file.write("    ")
                if not(self.parsedModel.ruleVariable[i] in self.parsedModel.speciesId):
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.ruleVariable[i]))+"]"
                    self.out_file.write(string)
                self.out_file.write("=")
    
                string = self.parsedModel.ruleFormula[i]
                for q in range(0,len(self.parsedModel.speciesId)):
                    pq = re.compile(self.parsedModel.speciesId[q])
                    string=pq.sub('y['+repr(q)+']' ,string)
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            string=pq.sub('tex2D(param_tex,'+repr(q)+',tid)',string)
    
                self.out_file.write(string)
                self.out_file.write(";\n")
                
        for i in range(0,len(self.parsedModel.listOfEvents)):
            self.out_file.write("    if( ")
            self.out_file.write(self.mathMLConditionParserCuda(self.parsedModel.eventCondition[i]))
            self.out_file.write("){\n")
            listOfAssignmentRules = self.parsedModel.listOfEvents[i].getListOfEventAssignments()
            for j in range(0, len(listOfAssignmentRules)):
                self.out_file.write("        ")
                #self.out_file.write("float ")
                if not(self.parsedModel.eventVariable[i][j] in self.parsedModel.speciesId):
                    self.out_file.write(self.parsedModel.eventVariable[i][j])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.eventVariable[i][j]))+"]"
                    self.out_file.write(string) 
                self.out_file.write("=")
                
                string = self.parsedModel.eventFormula[i][j]
                for q in range(0,len(self.parsedModel.speciesId)):
                    pq = re.compile(self.parsedModel.speciesId[q])
                    string=pq.sub('y['+repr(q)+']' ,string)
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            string=pq.sub('tex2D(param_tex,'+repr(q)+',tid)' ,string)
    
                self.out_file.write(string)
                self.out_file.write(";\n")
            self.out_file.write("}\n")
    
        self.out_file.write("\n")
    
        for i in range(0, len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isAssignment():
                self.out_file.write("    ")
                if not(self.parsedModel.ruleVariable[i] in self.parsedModel.speciesId):
                    self.out_file.write("float ")
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.ruleVariable[i]))+"]"
                    self.out_file.write(string)
                self.out_file.write("=")
     
                string = self.mathMLConditionParserCuda(self.parsedModel.ruleFormula[i])
                for q in range(0,len(self.parsedModel.speciesId)):
                    pq = re.compile(self.parsedModel.speciesId[q])
                    string=pq.sub("y["+repr(q)+"]" ,string)
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            x = "tex2D(param_tex,"+repr(q)+",tid)"
                            string=pq.sub(x,string)
                self.out_file.write(string)
                self.out_file.write(";\n")
        self.out_file.write("")
    
    
        #Write the derivatives
        for i in range(self.parsedModel.numSpecies-1,-1, -1):
            
            if (self.parsedModel.species[i].getConstant() == False and self.parsedModel.species[i].getBoundaryCondition() == False):
                self.out_file.write("    float d_y"+repr(i)+"= DT * (")
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write("(")
                
                reactionWritten = False
                for k in range(0,self.parsedModel.numReactions):
                    if(not self.parsedModel.stoichiometricMatrix[i][k]==0.0):
    
                        if(reactionWritten and self.parsedModel.stoichiometricMatrix[i][k]>0.0):
                            self.out_file.write("+")
                        reactionWritten = True
                        self.out_file.write(repr(self.parsedModel.stoichiometricMatrix[i][k]))
                        self.out_file.write("*(")
                        
                        #test if reaction has a positive sign
                        #if(reactionWritten):
                        #    if(stoichiometricMatrix[i][k]>0.0):
                        #        self.out_file.write("+")
                        #    else:
                        #        self.out_file.write("-")
                        #reactionWritten = True
                        
                        #test if reaction is 1.0; then omit multiplication term
                        #if(abs(stoichiometricMatrix[i][k]) == 1.0):
                        #    self.out_file.write("(")
                        #else:
                        #    self.out_file.write(repr(abs(stoichiometricMatrix[i][k])))
                        #    self.out_file.write("*(")
                            
                        string = self.parsedModel.kineticLaw[k]
                        for q in range(len(self.parsedModel.speciesId)-1,-1,-1):  
                            pq = re.compile(self.parsedModel.speciesId[q])
                            string=pq.sub('y['+repr(q)+']' ,string)
                        for q in range(0,len(self.parsedModel.parameterId)):
                            if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    pq = re.compile(self.parsedModel.parameterId[q])
                                    string=pq.sub('tex2D(param_tex,'+repr(q)+',tid)' ,string)
       
                        string=p.sub('',string)
                        
                        self.out_file.write(string)
                        self.out_file.write(")")
                        
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write(")/")
                    mySpeciesCompartment = self.parsedModel.species[i].getCompartment()
                    for j in range(0, len(self.parsedModel.listOfParameter)):
                        if (self.parsedModel.listOfParameter[j].getId() == mySpeciesCompartment):
                            if (not(self.parsedModel.parameterId[j] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[j] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    self.out_file.write("tex2D(param_tex,"+repr(j)+",tid)"+");")
                                    break
                                else:
                                    self.out_file.write(self.parsedModel.parameterId[j]+");")
                                    break
                else:
                    self.out_file.write(");")
                self.out_file.write("\n")
        
        self.out_file.write("\n")
    
        # check for columns of the stochiometry matrix with more than one entry
        randomVariables = ["*randNormal(rngRegs,sqrt(DT))"] * self.parsedModel.numReactions
        for k in range(0,self.parsedModel.numReactions):
            countEntries = 0
            for i in range(0,self.parsedModel.numSpecies):
                if(self.parsedModel.stoichiometricMatrix[i][k] != 0.0): countEntries += 1
            
            # define specific randomVariable
            if countEntries > 1:
                self.out_file.write("    float rand"+repr(k)+" = randNormal(rngRegs,sqrt(DT));\n")
                randomVariables[k] = "*rand" + repr(k)
        
        self.out_file.write("\n")
                
        #write noise terms
        for i in range(0,self.parsedModel.numSpecies):
            if (self.parsedModel.species[i].getConstant() == False and self.parsedModel.species[i].getBoundaryCondition() == False):
                self.out_file.write("    d_y"+repr(i)+" += (")
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write("(")
                
                reactionWritten = False
                for k in range(0,self.parsedModel.numReactions):
                    if(not self.parsedModel.stoichiometricMatrix[i][k]==0.0):
    
                        if(reactionWritten and self.parsedModel.stoichiometricMatrix[i][k]>0.0):
                            self.out_file.write("+")
                        reactionWritten = True
                        self.out_file.write(repr(self.parsedModel.stoichiometricMatrix[i][k]))
                        self.out_file.write("*sqrt(")
                        
                        #test if reaction has a positive sign
                        #if(reactionWritten):
                        #    if(stoichiometricMatrix[i][k]>0.0):
                        #        self.out_file.write("+")
                        #    else:
                        #        self.out_file.write("-")
                        #reactionWritten = True
                             
                        #test if reaction is 1.0; then omit multiplication term
                        #if(abs(stoichiometricMatrix[i][k]) == 1.0):
                        #    self.out_file.write("sqrtf(")
                        #else:
                        #    self.out_file.write(repr(abs(stoichiometricMatrix[i][k])))
                        #    self.out_file.write("*sqrtf(")
    
                        string = self.parsedModel.kineticLaw[k]
                        for q in range(0,len(self.parsedModel.speciesId)):
                            pq = re.compile(self.parsedModel.speciesId[q])
                            string=pq.sub('y['+repr(q)+']' ,string)
                        for q in range(0,len(self.parsedModel.parameterId)):
                            if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    pq = re.compile(self.parsedModel.parameterId[q])
                                    string=pq.sub('tex2D(param_tex,'+repr(q)+',tid)' ,string)
       
                        string=p.sub('',string)
                        self.out_file.write(string)
                        
                        #multiply random variable
                        self.out_file.write(")")
                        self.out_file.write(randomVariables[k])
                        #self.out_file.write("*randNormal(rngRegs,sqrt(DT))")
                        
                        
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write(")/")
                    mySpeciesCompartment = self.parsedModel.species[i].getCompartment()
                    for j in range(0, len(self.parsedModel.listOfParameter)):
                        if (self.parsedModel.listOfParameter[j].getId() == mySpeciesCompartment):
                            if (not(self.parsedModel.parameterId[j] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[j] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    self.out_file.write("tex2D(param_tex,"+repr(j)+",tid)"+")")
                                    break
                                else:
                                    self.out_file.write(self.parsedModel.parameterId[j]+")")
                                    break
                else:
                    self.out_file.write(")")
                self.out_file.write(";\n")
        
        self.out_file.write("\n")
        #add terms
        for i in range(0,self.parsedModel.numSpecies):
            if (self.parsedModel.species[i].getConstant() == False and self.parsedModel.species[i].getBoundaryCondition() == False ):
                self.out_file.write("    y["+repr(i)+"] += d_y"+repr(i)+";\n")
            
        self.out_file.write("}\n")
    
        
    ################# same file
    
    
        p=re.compile('\s')
        #The user-defined functions used in the model must be written in the file
        self.out_file.write("//Code for shared memory\n")
    
        numEvents = len(self.parsedModel.listOfEvents)
        numRules = len(self.parsedModel.listOfRules)
        num = numEvents+numRules
        if num>0:
            self.out_file.write("#define leq(a,b) a<=b\n")
            self.out_file.write("#define neq(a,b) a!=b\n")
            self.out_file.write("#define geq(a,b) a>=b\n")
            self.out_file.write("#define lt(a,b) a<b\n")
            self.out_file.write("#define gt(a,b) a>b\n")
            self.out_file.write("#define eq(a,b) a==b\n")
            self.out_file.write("#define and_(a,b) a&&b\n")
            self.out_file.write("#define or_(a,b) a||b\n")
    
        for i in range(0,len(self.parsedModel.listOfFunctions)):
            self.out_file.write("__device__ float "+self.parsedModel.listOfFunctions[i].getId()+"(")
            for j in range(0, self.parsedModel.listOfFunctions[i].getNumArguments()):
                self.out_file.write("float "+self.parsedModel.functionArgument[i][j])
                if(j<( self.parsedModel.listOfFunctions[i].getNumArguments()-1)):
                    self.out_file.write(",")
            self.out_file.write("){\n    return ")
            self.out_file.write(self.parsedModel.functionBody[i])
            self.out_file.write(";\n}\n")
            self.out_file.write("\n")
    
    
    
    
    
    
        self.out_file.write("\n")
        self.out_file.write("__device__ void step(float *parameter, float *y, float t, unsigned *rngRegs){\n")
    
        numSpecies = len(self.parsedModel.species)
    
        #write rules and events
        for i in range(0,len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isRate() == True:
                self.out_file.write("    ")
                if not(self.parsedModel.ruleVariable[i] in self.parsedModel.speciesId):
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(ruleVariable[i]))+"]"
                    self.out_file.write(string)
                self.out_file.write("=")
    
                string = self.parsedModel.ruleFormula[i]
                for q in range(0,len(self.parsedModel.speciesId)):
                    pq = re.compile(self.parsedModel.speciesId[q])
                    string=pq.sub('y['+repr(q)+']' ,string)
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(EventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            string=pq.sub('parameter['+repr(q)+']' ,string)
    
                self.out_file.write(string)
                self.out_file.write(";\n")
                
        for i in range(0,len(self.parsedModel.listOfEvents)):
            self.out_file.write("    if( ")
            self.out_file.write(self.mathMLConditionParserCuda(self.parsedModel.eventCondition[i]))
            self.out_file.write("){\n")
            listOfAssignmentRules = self.parsedModel.listOfEvents[i].getListOfEventAssignments()
            for j in range(0, len(listOfAssignmentRules)):
                self.out_file.write("        ")
                #self.out_file.write("float ")
                if not(self.parsedModel.eventVariable[i][j] in self.parsedModel.speciesId):
                    self.out_file.write(self.parsedModel.eventVariable[i][j])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.eventVariable[i][j]))+"]"
                    self.out_file.write(string) 
                self.out_file.write("=")
                
                string = self.parsedModel.eventFormula[i][j]
                for q in range(0,len(self.parsedModel.speciesId)):
                    pq = re.compile(self.parsedModel.speciesId[q])
                    string=pq.sub('y['+repr(q)+']' ,string)
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            string=pq.sub('parameter['+repr(q)+']' ,string)
    
                self.out_file.write(string)
                self.out_file.write(";\n")
            self.out_file.write("}\n")
    
        self.out_file.write("\n")
    
        for i in range(0, len(self.parsedModel.listOfRules)):
            if self.parsedModel.listOfRules[i].isAssignment():
                self.out_file.write("    ")
                if not(self.parsedModel.ruleVariable[i] in self.parsedModel.speciesId):
                    self.out_file.write("float ")
                    self.out_file.write(self.parsedModel.ruleVariable[i])
                else:
                    string = "y["+repr(self.parsedModel.speciesId.index(self.parsedModel.ruleVariable[i]))+"]"
                    self.out_file.write(string)
                self.out_file.write("=")
     
                string = self.mathMLConditionParserCuda(self.parsedModel.ruleFormula[i])
                for q in range(0,len(self.parsedModel.speciesId)):
                    pq = re.compile(self.parsedModel.speciesId[q])
                    string=pq.sub("y["+repr(q)+"]" ,string)
                for q in range(0,len(self.parsedModel.parameterId)):
                    if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                        flag = False
                        for r in range(0,len(self.parsedModel.eventVariable)):
                            if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                flag = True
                        if flag==False:
                            pq = re.compile(self.parsedModel.parameterId[q])
                            x = "parameter["+repr(q)+"]"
                            string=pq.sub(x,string)
                self.out_file.write(string)
                self.out_file.write(";\n")
        #self.out_file.write("\n\n")
    
    
        #Write the derivatives
        for i in range(0,self.parsedModel.numSpecies):
            if (self.parsedModel.species[i].getConstant() == False and self.parsedModel.species[i].getBoundaryCondition() == False):
                self.out_file.write("    float d_y"+repr(i)+"= DT * (")
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write("(")
                
                reactionWritten = False
                for k in range(0,self.parsedModel.numReactions):
                    if(not self.parsedModel.stoichiometricMatrix[i][k]==0.0):
    
                        if(reactionWritten and self.parsedModel.stoichiometricMatrix[i][k]>0.0):
                            self.out_file.write("+")
                        reactionWritten = True
                        self.out_file.write(repr(self.parsedModel.stoichiometricMatrix[i][k]))
                        self.out_file.write("*(")
                        
                        #test if reaction has a positive sign
                        #if(reactionWritten):
                        #    if(stoichiometricMatrix[i][k]>0.0):
                        #        self.out_file.write("+")
                        #    else:
                        #        self.out_file.write("-")
                        #reactionWritten = True
                        
                        #test if reaction is 1.0; then omit multiplication term
                        #if(abs(stoichiometricMatrix[i][k]) == 1.0):
                        #    self.out_file.write("(")
                        #else:
                        #    self.out_file.write(repr(abs(stoichiometricMatrix[i][k])))
                        #    self.out_file.write("*(")
    
                        string = self.parsedModel.kineticLaw[k]
                        for q in range(0,len(self.parsedModel.speciesId)):
                            pq = re.compile(self.parsedModel.speciesId[q])
                            string=pq.sub('y['+repr(q)+']' ,string)
                        for q in range(0,len(self.parsedModel.parameterId)):
                            if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    pq = re.compile(self.parsedModel.parameterId[q])
                                    string=pq.sub('parameter['+repr(q)+']' ,string)
       
                        string=p.sub('',string)
                        
                        self.out_file.write(string)
                        self.out_file.write(")")
                        
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write(")/")
                    mySpeciesCompartment = self.parsedModel.species[i].getCompartment()
                    for j in range(0, len(self.parsedModel.listOfParameter)):
                        if (self.parsedModel.listOfParameter[j].getId() == mySpeciesCompartment):
                            if (not(self.parsedModel.parameterId[j] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[j] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    self.out_file.write("parameter["+repr(j)+"]"+");")
                                    break
                                else:
                                    self.out_file.write(self.parsedModel.parameterId[j]+");")                                
                                    break
                else:
                    self.out_file.write(");")
                self.out_file.write("\n")
    
        self.out_file.write("\n")
    
        # check for columns of the stochiometry matrix with more than one entry
        randomVariables = ["*randNormal(rngRegs,sqrt(DT))"] * self.parsedModel.numReactions
        for k in range(0,self.parsedModel.numReactions):
            countEntries = 0
            for i in range(0,numSpecies):
                if(self.parsedModel.stoichiometricMatrix[i][k] != 0.0): countEntries += 1
            
            # define specific randomVariable
            if countEntries > 1:
                self.out_file.write("    float rand"+repr(k)+" = randNormal(rngRegs,sqrt(DT));\n")
                randomVariables[k] = "*rand" + repr(k)
        
        self.out_file.write("\n")
    
        #write noise terms
        for i in range(0,self.parsedModel.numSpecies):
            if (self.parsedModel.species[i].getConstant() == False and self.parsedModel.species[i].getBoundaryCondition() == False):
                self.out_file.write("    d_y"+repr(i)+"+= (")
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write("(")
                    
                reactionWritten = False
                for k in range(0,self.parsedModel.numReactions):
                    if(not self.parsedModel.stoichiometricMatrix[i][k]==0.0):
    
                        if(reactionWritten and self.parsedModel.stoichiometricMatrix[i][k]>0.0):
                            self.out_file.write("+")
                        reactionWritten = True
                        self.out_file.write(repr(self.parsedModel.stoichiometricMatrix[i][k]))
                        self.out_file.write("*sqrt(")
                        
                        #test if reaction has a positive sign
                        #if(reactionWritten):
                        #    if(stoichiometricMatrix[i][k]>0.0):
                        #        self.out_file.write("+")
                        #    else:
                        #        self.out_file.write("-")
                        #reactionWritten = True
                             
                        #test if reaction is 1.0; then omit multiplication term
                        #if(abs(stoichiometricMatrix[i][k]) == 1.0):
                        #    self.out_file.write("sqrtf(")
                        #else:
                        #    self.out_file.write(repr(abs(stoichiometricMatrix[i][k])))
                        #    self.out_file.write("*sqrtf(")
    
                        string = self.parsedModel.kineticLaw[k]
                        for q in range(0,len(self.parsedModel.speciesId)):
                            pq = re.compile(self.parsedModel.speciesId[q])
                            string=pq.sub('y['+repr(q)+']' ,string)
                        for q in range(0,len(self.parsedModel.parameterId)):
                            if (not(self.parsedModel.parameterId[q] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[q] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    pq = re.compile(self.parsedModel.parameterId[q])
                                    string=pq.sub('parameter['+repr(q)+']' ,string)
       
                        string=p.sub('',string)
                        self.out_file.write(string)
                        
                        #multiply random variable
                        self.out_file.write(")")
                        self.out_file.write(randomVariables[k])
                        #self.out_file.write("*randNormal(rngRegs,sqrt(DT))")
                        
                if (self.parsedModel.species[i].isSetCompartment() == True):
                    self.out_file.write(")/")
                    mySpeciesCompartment = self.parsedModel.species[i].getCompartment()
                    for j in range(0, len(self.parsedModel.listOfParameter)):
                        if (self.parsedModel.listOfParameter[j].getId() == mySpeciesCompartment):
                            if (not(self.parsedModel.parameterId[j] in self.parsedModel.ruleVariable)):
                                flag = False
                                for r in range(0,len(self.parsedModel.eventVariable)):
                                    if (self.parsedModel.parameterId[j] in self.parsedModel.eventVariable[r]):
                                        flag = True
                                if flag==False:
                                    self.out_file.write("parameter["+repr(j)+"]"+")")
                                    break
                                else:
                                    self.out_file.write(self.parsedModel.parameterId[j]+")")                                
                                    break
                else:
                    self.out_file.write(")")
    
                self.out_file.write(";\n")
    
        self.out_file.write("\n")
        #add terms
        for i in range(0,self.parsedModel.numSpecies):
            if (self.parsedModel.species[i].getConstant() == False and self.parsedModel.species[i].getBoundaryCondition() == False):
                self.out_file.write("    y["+repr(i)+"] += d_y"+repr(i)+";\n")
            
        self.out_file.write("}\n")
