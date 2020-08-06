#include "ParameterSet.h"

#include <iostream>

std::size_t ParameterSet::addParameter(std::string name){
	if(nameMap.count(name))
		throw std::logic_error("A parameter with the name '"+name+
							   "' already exists in the is ParameterSet");
	std::size_t idx=names.size();
	names.push_back(name);
	parameterValues.push_back(0);
	lowerBounds.push_back(-std::numeric_limits<double>::infinity());
	upperBounds.push_back(std::numeric_limits<double>::infinity());
	//properties.push_back(PropertyCollection());
	nameMap.emplace(name,idx);
	return(idx);
}
	
std::size_t ParameterSet::getParameterIndex(std::string name) const{
	auto it=nameMap.find(name);
	if(it==nameMap.end())
		throw std::logic_error("No parameter with the name '"+name+
							   "' exists in the is ParameterSet");
	return(it->second);
}
	
const std::string& ParameterSet::getParameterName(std::size_t index) const{
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	return(names[index]);
}
	
std::size_t ParameterSet::getFreeParameterIndex(std::size_t freeIndex) const{
	if(freeIndex>=numberOfFreeParameters())
		throw std::logic_error("Index "+std::to_string(freeIndex)+" is out of range");
	std::size_t index=0;
	for(std::size_t j=0, k=0; ; index++){
		if(j<fixedParams.size() && fixedParams[j]==index)
			j++;
		else
			k++;
		if(k==freeIndex+1)
			break;
	}
	return(index);
}
	
void ParameterSet::setParameterValue(std::size_t index, double value){
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	parameterValues[index]=value;
}
	
void ParameterSet::setParameterValue(const std::string& name, double value){
	parameterValues[getParameterIndex(name)]=value;
}

double ParameterSet::getParameterValue(std::size_t index) const{
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	return(parameterValues[index]);
}

double ParameterSet::getParameterValue(const std::string& name) const{
	return(parameterValues[getParameterIndex(name)]);
}
	
void ParameterSet::setParameterLowerLimit(std::size_t index, double limit){
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	lowerBounds[index]=limit;
}
	
void ParameterSet::setParameterLowerLimit(const std::string& name, double limit){
	setParameterLowerLimit(getParameterIndex(name),limit);
}
	
double ParameterSet::getParameterLowerLimit(std::size_t index) const{
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	return(lowerBounds[index]);
}
	
double ParameterSet::getParameterLowerLimit(const std::string& name) const{
	return(getParameterLowerLimit(getParameterIndex(name)));
}
	
void ParameterSet::setParameterUpperLimit(std::size_t index, double limit){
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	upperBounds[index]=limit;
}
	
void ParameterSet::setParameterUpperLimit(const std::string& name, double limit){
	setParameterUpperLimit(getParameterIndex(name),limit);
}
	
double ParameterSet::getParameterUpperLimit(std::size_t index) const{
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	return(upperBounds[index]);
}
	
double ParameterSet::getParameterUpperLimit(const std::string& name) const{
	return(getParameterUpperLimit(getParameterIndex(name)));
}
	
void ParameterSet::fixParameter(size_t index){
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	if(std::find(fixedParams.begin(),fixedParams.end(),index)==fixedParams.end()){
		if(std::binary_search(fixedParams.begin(),fixedParams.end(),index))
			return;
		fixedParams.push_back(index);
		std::sort(fixedParams.begin(),fixedParams.end());
	}
}
	
void ParameterSet::fixParameter(const std::string& name){
	fixParameter(getParameterIndex(name));
}
	
void ParameterSet::freeParameter(size_t index){
	auto it=std::find(fixedParams.begin(),fixedParams.end(),index);
	if(it!=fixedParams.end())
		fixedParams.erase(it);
}
	
void ParameterSet::freeParameter(const std::string& name){
	freeParameter(getParameterIndex(name));
}
	
bool ParameterSet::isFixed(std::size_t idx) const{
	if(idx>=names.size())
		throw std::logic_error("Index "+std::to_string(idx)+" is out of range");
	return(std::find(fixedParams.begin(),fixedParams.end(),idx)!=fixedParams.end());
}
	
bool ParameterSet::isFixed(const std::string& name) const{
	return(isFixed(getParameterIndex(name)));
}

std::size_t ParameterSet::numberOfParameters() const{
	return(names.size());
}

std::size_t ParameterSet::numberOfFreeParameters() const{
	return(numberOfParameters()-numberOfFixedParameters());
}

std::size_t ParameterSet::numberOfFixedParameters() const{
	return(fixedParams.size());
}
	
const std::vector<double>& ParameterSet::getParameterValues() const{
	return(parameterValues);
}
	
std::vector<double> ParameterSet::getFreeParameterValues() const{
	std::vector<double> values(numberOfFreeParameters());
	externalToInternal(parameterValues,values);
	return(values);
}
	
std::vector<double> ParameterSet::getFreeParameterLowerBounds() const{
	std::vector<double> lbounds(numberOfFreeParameters());
	externalToInternal(lowerBounds,lbounds);
	return(lbounds);
}
	
std::vector<double> ParameterSet::getFreeParameterUpperBounds() const{
	std::vector<double> ubounds(numberOfFreeParameters());
	externalToInternal(upperBounds,ubounds);
	return(ubounds);
}
	
void ParameterSet::setParameterValues(const std::vector<double>& values){
	if(values.size()!=numberOfParameters())
		throw std::logic_error("Incorrect number of parameter values: "
							   +std::to_string(values.size())+", should be "
							   +std::to_string(numberOfParameters()));
	parameterValues=values;
}
	
void ParameterSet::setFreeParameterValues(const std::vector<double>& values){
	if(values.size()!=numberOfFreeParameters())
		throw std::logic_error("Incorrect number of free parameter values: "
							   +std::to_string(values.size())+", should be "
							   +std::to_string(numberOfFreeParameters()));
	internalToExternal(values,parameterValues);
}
	
void ParameterSet::insertFreeParameters(const std::vector<double>& freeValues, std::vector<double>& allValues) const{
	if(freeValues.size()!=numberOfFreeParameters())
		throw std::logic_error("Incorrect number of free parameter values: "
							   +std::to_string(freeValues.size())+", should be "
							   +std::to_string(numberOfFreeParameters()));
	if(allValues.size()!=numberOfParameters())
		throw std::logic_error("Incorrect number of parameter values: "
							   +std::to_string(allValues.size())+", should be "
							   +std::to_string(numberOfParameters()));
	internalToExternal(freeValues,allValues);
}

void ParameterSet::extractFreeParameters(const std::vector<double>& allValues, std::vector<double>& freeValues) const{
	if(freeValues.size()!=numberOfFreeParameters())
		throw std::logic_error("Incorrect number of free parameter values: "
							   +std::to_string(freeValues.size())+", should be "
							   +std::to_string(numberOfFreeParameters()));
	if(allValues.size()!=numberOfParameters())
		throw std::logic_error("Incorrect number of parameter values: "
							   +std::to_string(allValues.size())+", should be "
							   +std::to_string(numberOfParameters()));
	externalToInternal(allValues,freeValues);
}
	
std::vector<double> ParameterSet::collateValues(const std::map<std::string,double>& values){
	std::vector<double> result=getParameterValues();
	for(const auto& param : values)
		result[getParameterIndex(param.first)]=param.second;
	return(result);
}
	
bool ParameterSet::inBounds(std::size_t index, double value) const{
	if(index>=names.size())
		throw std::logic_error("Index "+std::to_string(index)+" is out of range");
	return(value>=lowerBounds[index] && value<=upperBounds[index]);
}
	
bool ParameterSet::inBounds(const std::string name, double value) const{
	return(inBounds(getParameterIndex(name),value));
}
	
bool ParameterSet::inBounds(const std::vector<double>& values) const{
	if(values.size()!=numberOfParameters())
		throw std::logic_error("Incorrect number of parameter values: "
							   +std::to_string(values.size())+", should be "
							   +std::to_string(numberOfParameters()));
	for(std::size_t i=0; i<numberOfParameters(); i++){
		if(!inBounds(i,values[i]))
			return(false);
	}
	return(true);
}

