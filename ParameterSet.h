#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H

#include <algorithm>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
	
class ParameterSet{
public:
	///Add a new parameter with the given name.
	///The parameter's default value will be zero, and it will be free and 
	///unbounded.
	///\param name the name of the parameter to create
	///\return the index of the newly created parameter
	///\throws std::logic_error if a parameter with the given name already exists
	std::size_t addParameter(std::string name);
	///Obtain the index of the parameter with a given name
	///\param name the name of the parameter to look up
	///\throws std::logic_error if no parameter with the given name exists
	std::size_t getParameterIndex(std::string name) const;
	///Obtain the name of the parameter with a given index
	///\param index the index of the parameter to look fetch
	///\throws std::logic_error if the given index is out of range
	const std::string& getParameterName(std::size_t index) const;
	///Obtain the overall index of the parameter at a given index among the 
	///currently free parameters
	///\param freeIndex index within the set of free parameters
	///\throws std::logic_error if the given index is out of range
	std::size_t getFreeParameterIndex(std::size_t freeIndex) const;
	
	///Set the initial/current value of one parameter
	///\param index the index of the parameter
	///\param value the new value for the parameter
	///\throws std::logic_error if the given index is out of range
	void setParameterValue(std::size_t index, double value);
	///Set the initial/current value of one parameter
	///\param name the name of the parameter to look up
	///\param value the new value for the parameter
	///\throws std::logic_error if no parameter with the given name exists
	void setParameterValue(const std::string& name, double value);
	
	///Get the initial/current value of one parameter
	///\param index the index of the parameter
	///\throws std::logic_error if the given index is out of range
	double getParameterValue(std::size_t index) const;
	///Get the initial/current value of one parameter
	///\param name the name of the parameter
	///\throws std::logic_error if no parameter with the given name exists
	double getParameterValue(const std::string& name) const;
	
	///Set the lower limit allowed for one parameter
	///\param index the index of the parameter
	///\param limit the limit on the parameter. If this is -infinity any
	///             existing limit is removed.
	///\throws std::logic_error if the given index is out of range
	void setParameterLowerLimit(std::size_t index, double limit);
	///Set the lower limit allowed for one parameter
	///\param name the name of the parameter
	///\param limit the limit on the parameter. If this is -infinity any
	///             existing limit is removed.
	///\throws std::logic_error if no parameter with the given name exists
	void setParameterLowerLimit(const std::string& name, double limit);
	
	///Get the lower limit allowed for one parameter
	///\param index the index of the parameter
	///\throws std::logic_error if the given index is out of range
	double getParameterLowerLimit(std::size_t index) const;
	///Get the lower limit allowed for one parameter
	///\param name the name of the parameter
	///\throws std::logic_error if no parameter with the given name exists
	double getParameterLowerLimit(const std::string& name) const;
	
	///Set the upper limit allowed for one parameter
	///\param index the index of the parameter
	///\param limit the limit on the parameter. If this is infinity any
	///             existing limit is removed.
	///\throws std::logic_error if the given index is out of range
	void setParameterUpperLimit(std::size_t index, double limit);
	///Set the lower limit allowed for one parameter
	///\param name the name of the parameter
	///\param limit the limit on the parameter. If this is infinity any
	///             existing limit is removed.
	///\throws std::logic_error if no parameter with the given name exists
	void setParameterUpperLimit(const std::string& name, double limit);
	
	///Get the upper limit allowed for one parameter
	///\param index the index of the parameter
	///\throws std::logic_error if the given index is out of range
	double getParameterUpperLimit(std::size_t index) const;
	///Get the lower limit allowed for one parameter
	///\param name the name of the parameter
	///\throws std::logic_error if no parameter with the given name exists
	double getParameterUpperLimit(const std::string& name) const;
	
	///Hold a parameter fixed.
	///This ignores and leaves unchanged any bounds which may be set for the 
	///parameter.
	///\param index the index of the parameter
	///\throws std::logic_error if the given index is out of range
	void fixParameter(std::size_t index);
	///Hold a parameter fixed.
	///This ignores and leaves unchanged any bounds which may be set for the 
	///parameter.
	///\param name the name of the parameter
	///\throws std::logic_error if no parameter with the given name exists
	void fixParameter(const std::string& name);
	///Free a fixed parameter
	///This ignores and leaves unchanged any bounds which may be set for the 
	///parameter. If the parameter is already free, this function has no effect.
	///\param index the index of the parameter
	///\throws std::logic_error if the given index is out of range
	void freeParameter(std::size_t index);
	///Free a fixed parameter
	///This ignores and leaves unchanged any bounds which may be set for the 
	///parameter. If the parameter is already free, this function has no effect.
	///\param name the name of the parameter
	///\throws std::logic_error if no parameter with the given name exists
	void freeParameter(const std::string& name);
	
	///Determine whether a given parameter is currently held fixed
	///\param index the index of the parameter
	///\throws std::logic_error if the given index is out of range
	bool isFixed(std::size_t index) const;
	///Determine whether a given parameter is currently held fixed
	///\param name the name of the parameter
	///\throws std::logic_error if no parameter with the given name exists
	bool isFixed(const std::string& name) const;
	
	///Obtain the total number of parameters
	std::size_t numberOfParameters() const;
	///Obtain the number of parameters which are currently free
	std::size_t numberOfFreeParameters() const;
	///Obtain the number of parameters which are currently fixed
	std::size_t numberOfFixedParameters() const;
	
	///Obtain the current values of all parameters, in order.
	const std::vector<double>& getParameterValues() const;
	///Obtain the values of all currently free parameters, in order. 
	std::vector<double> getFreeParameterValues() const;
	///Obtain the lower bounds for all currently free parameters, in order. 
	std::vector<double> getFreeParameterLowerBounds() const;
	///Obtain the upper bounds for all currently free parameters, in order. 
	std::vector<double> getFreeParameterUpperBounds() const;
	///Set the values of all parameters
	///\param value the ordered collection of new values to use
	void setParameterValues(const std::vector<double>& values);
	///Set the values of all currently free parameters
	///\param value the ordered collection of new values to use
	void setFreeParameterValues(const std::vector<double>& values);
	
	///Copy a set of new values for the free parameters into a complete list of
	///all parameter values
	///\param freeValues the ordered set of new values for the free parameters
	///\param allValues the full ordered set of parameter values into which the 
	///                 new values of the free parameters will be inserted
	void insertFreeParameters(const std::vector<double>& freeValues, std::vector<double>& allValues) const;
	
	///Copy the values associated with free parameters out of a set of values for
	///all parameters
	///\param allValues the full collection of parameter values from which to 
	///                 extract those corresponding to free parameters
	///\param freeValues the location where the free parameter values will be stored
	void extractFreeParameters(const std::vector<double>& allValues, std::vector<double>& freeValues) const;
	
	///Extract one value from an ordered collection of values
	///\param name the name of the parameter to fetch
	///\param values the full collection of values
	///\pre values.size()==numberOfParameters()
	///\throws std::logic_error if no parameter with the given name exists
	template<typename Container, typename T=typename Container::value_type>
	const T& extractParameter(const std::string name, const Container& values){
		return(values[getParameterIndex(name)]);
	}
	
	///Sort a map of parameter names to values into a vector. 
	///Parameters which have no value specified in the map will use the default
	///value previously stored in the ParameterSet.
	///\param values the non-default parameter value to use
	///\throws std::logic_error if any key in the map does not match the name of
	///                         some parameter
	std::vector<double> collateValues(const std::map<std::string,double>& values);
	
	///Test parameter values is within allowed bounds. 
	///NaN values will always fail.
	///\param index the index of the parameter
	///\param value the value to be tested against the parameter's bounds
	///\throws std::logic_error if the given index is out of range
	bool inBounds(std::size_t index, double value) const;
	
	///Test parameter values is within allowed bounds. 
	///NaN values will always fail.
	///\param index the index of the parameter
	///\param value the value to be tested against the parameter's bounds
	///\throws std::logic_error if no parameter with the given name exists
	bool inBounds(const std::string name, double value) const;
	
	///Test whether a set of parameter values are all within allowed bounds.
	///NaN values will always fail.
	///\param values the values to be tested against the parameters' bounds
	///\throws std::logic_error if the number of values given does not match the 
	///                         number of currently defined parameters
	bool inBounds(const std::vector<double>& values) const;
	
	///Associate a user-defined property with a parameter
	///\param index the index of the parameter
	///\param propName the name to use for the property
	///\param data the value for the property
	///\throws std::logic_error if the given index is out of range
	/*template<typename T>
	void setParameterProperty(std::size_t index, std::string propName, T data){
		if(index>=names.size())
			throw std::logic_error("Index "+std::to_string(index)+" is out of range");
		properties[index].setProperty(propName,boost::any(data));
	}*/
	///Associate a user-defined property with a parameter
	///\param name the name of the parameter
	///\param propName the name to use for the property
	///\param data the value for the property
	///\throws std::logic_error if no parameter with the given name exists
	/*template<typename T>
	void setParameterProperty(const std::string name, std::string propName, T data){
		properties[getParameterIndex(name)].setProperty(propName,boost::any(data));
	}*/
	///Retrieve a user-defined property associated with a parameter
	///\param index the index of the parameter
	///\param propName the name to use for the property
	///\throws std::logic_error if the given index is out of range, 
	///        std::logic_error if the property does not exist,
	///        boost::bad_any_cast if the proerty has a different type from T
	/*template<typename T>
	T getParameterProperty(std::size_t index, std::string propName) const{
		if(index>=names.size())
			throw std::logic_error("Index "+std::to_string(index)+" is out of range");
		return(boost::any_cast<T>(properties[index].getProperty(propName)));
	}*/
	///Retrieve a user-defined property associated with a parameter
	///\param name the name of the parameter
	///\param propName the name to use for the property
	///\throws std::logic_error if no parameter with the given name exists, 
	///        std::logic_error if the property does not exist,
	///        boost::bad_any_cast if the proerty has a different type from T
	template<typename T>
	T getParameterProperty(const std::string name, std::string propName) const{
		return(getParameterProperty<T>(getParameterIndex(name), propName));
	}
	///Check whether a parameter has a given user defined property associated with it
	///\param index the index of the parameter
	///\param propName the name to use for the property
	///\throws std::logic_error if the given index is out of range
	/*bool parameterHasProperty(std::size_t index, std::string propName) const{
		if(index>=names.size())
			throw std::logic_error("Index "+std::to_string(index)+" is out of range");
		return(properties[index].hasProperty(propName));
	}*/
	///Check whether a parameter has a given user defined property associated with it
	///\param name the name of the parameter
	///\param propName the name to use for the property
	///\throws std::logic_error if no parameter with the given name exists
	/*bool parameterHasProperty(const std::string name, std::string propName) const{
		return(parameterHasProperty(getParameterIndex(name), propName));
	}*/
	
private:
	std::vector<std::string> names;
	std::vector<double> parameterValues;
	std::vector<double> lowerBounds;
	std::vector<double> upperBounds;
	std::vector<size_t> fixedParams;
	
	/*struct PropertyCollection{
		std::unordered_map<std::string,boost::any> properties;
		
		void setProperty(const std::string& name, boost::any value){
			properties[name]=value;
		}
		bool hasProperty(const std::string& name) const{
			return(properties.count(name));
		}
		boost::any getProperty(const std::string& name) const{
			auto it=properties.find(name);
			if(it==properties.end())
				throw std::logic_error("No property named "+name);
			return(it->second);
		}
	};
	
	std::vector<PropertyCollection> properties;*/
	
	std::unordered_map<std::string,std::size_t> nameMap;
	
	///Copy values from a vector with all parameters to a vector including only 
	///free parameters
	///\param externalVector the evaulation position including all parameters
	///\param internalVector the evaulation position including only free parameters
	///\pre internalVector.size()==externalVector.size()-fixedParams.size()
	template<typename T>
	void externalToInternal(const std::vector<T>& externalVector,
							std::vector<T>& internalVector) const{
		//i is the external index
		//j is the internal index
		//k is the index in fixedParams
		const std::size_t nFixed=fixedParams.size();
		for(std::size_t i=0, j=0, k=0; i<numberOfParameters(); i++){
			if(k<nFixed && fixedParams[k]==i){
				k++;
				continue;
			}
			internalVector[j]=externalVector[i];
			j++;
		}
	}
	
	///Copy values from a vector with only free parameters to a vector including 
	///all parameters, without disturbing the
	///existing values in the full vector corresponding to fixed parameters
	///\param internalVector the evaulation position including only free parameters
	///\param externalVector the evaulation position including all parameters
	///\pre internalVector.size()==externalVector.size()-fixedParams.size()
	template<typename T>
	void internalToExternal(const std::vector<T>& internalVector,
							std::vector<T>& externalVector) const{
		//i is the external index
		//j is the internal index
		//k is the index in fixedParams
		const std::size_t nFixed=fixedParams.size();
		for(std::size_t i=0, j=0, k=0; i<numberOfParameters(); i++){
			if(k<nFixed && fixedParams[k]==i){
				k++;
				continue;
			}
			externalVector[i]=internalVector[j];
			j++;
		}
	}
};
	
#endif