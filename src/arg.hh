/**
   @file arg.hh
   @brief Template header file for command-line processing
*/
#pragma once
#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <cassert>
namespace arg {
	/// Base class for proxy to value of an Option
	class Value
	{
	public:
		virtual ~Value() {}
		virtual void set(std::string const & str) = 0; ///<convert the str to value and put it in storage
		virtual std::string to_str() const = 0; ///<convert the value to a string
		virtual std::string get_type() const = 0; ///<type name of the value
	};
	/// Value with a type
	template<typename T>
	class TValue :
		public Value
	{
		T & v;
	public:
		TValue(T & v) : v(v) {} ///<Construct a typed value
		void set(std::string const & str) override {(void)str;}
		std::string to_str() const override {return std::to_string(v);}
		std::string get_type() const override {return typeid(T).name();}
	};
	///@cond
	/// Specializations for typed values
	template<> void TValue<int>::set(std::string const & str) {v = std::stoi(str);}
	template<> void TValue<double>::set(std::string const & str) {v = std::stod(str);}
	template<> void TValue<std::string>::set(std::string const & str) {v=str;}
	template<> std::string TValue<std::string>::to_str() const {return v;}
	template<> std::string TValue<int>::get_type() const {return "int";}
	template<> std::string TValue<unsigned>::get_type() const {return "uint";}
	template<> std::string TValue<double>::get_type() const {return "real";}
	template<> std::string TValue<std::string>::get_type() const {return "string";}
	///@endcond

	/// Value with a storage
	class ValueStore
	{
	protected:
		std::shared_ptr<Value> store_ptr; ///<Value storage
		std::string _desc; ///<Description storage
	public:
		std::string type_label() const ///<Get a type label for the value
		{
			if (store_ptr) {
				auto l = store_ptr->get_type();
				for (auto &c:l) c = toupper(c);
				return l;
			}
			return "?";
		}
		std::string get_desc() const {return _desc;} ///<Get description of the value
	};

	/// Command-line parameter
	class Param :
		public ValueStore
	{
		std::string _name; // Parameter name
	public:
		/// Set value of the parameter to string
		/// @param str String representing the value
		void set_value(std::string const & str) {store_ptr->set(str);}
		std::string get_str() const {return store_ptr->to_str();} ///<Get value of the parameter as a string
		std::string get_name() const {return _name;} ///<Get name of the parameter
		std::string get_label() const {return _name.length()?_name:type_label();} ///<Get label of the parameter
		///\name Parameter modifiers
		/// Chaining modifiers for a parameter
		///@{
		Param &store(std::shared_ptr<Value> p) {store_ptr = p;return *this;} ///<Where to store
		Param &name(std::string const &n) {_name = n;return *this;} ///<Name the parameter
		Param &desc(std::string const &s) {_desc = s;return *this;} ///<Describe the parameter
		///@}
		friend class Parser;
	};

	/// A command-line option
	class Option :
		public ValueStore
	{
		std::string _id; // Option id string
		std::string _def; // default value
		bool _present; // present?
	public:
		/// Construct option with given id string
		Option(std::string const &n) : _id(n),vr(VR_YES) {}
		/// Set value of the option to string
		/// @param str String representing the option value
		void set_value(std::string const & str) {store_ptr->set(str);}
		/// Set value of the option to default
		void set_value() {store_ptr->set(_def);}
		/// Get id of the option
		std::string const &get_id() const {return _id;}
		/// Value requirement types
		enum ValueReq {VR_YES,VR_NO,VR_MAYBE};
		int vr; ///<Value requirement
		///\name Option modifiers
		/// Chaining modifies for an option
		///@{
		Option &store(std::shared_ptr<Value> p) {store_ptr = p;return *this;} ///<Where to store
		Option &id(std::string const &i) {_id = i;return *this;} ///<Set id
		Option &fill(std::string const &s) {_def = s;vr = VR_MAYBE;return *this;} ///<Optional value
		Option &fix(std::string const &s) {_def = s;vr = VR_NO;return *this;} ///<No value
		Option &desc(std::string const &s) {_desc = s;return *this;} ///<Description
		///@}
		friend class Parser;
	};

	/// The command-line parser
	class Parser
	{
		std::vector<std::shared_ptr<Param>> prm_list;
		unsigned prm_required; /// number of required parameters
		std::vector<std::shared_ptr<Option>> opt_list;
		std::string last_cmd;
	public:
		Parser() : prm_required(0) {}
		~Parser() {}
		/// Add a parameter
		template<typename T>
		Param &add_prm(T &v)
		{
			auto p = std::make_shared<Param>();
			p->store(std::make_shared<TValue<T>>(v));
			prm_list.push_back(p);
			return *p;
		}
		/// Paramters added so far are required
		void required() {prm_required = prm_list.size();} // set number of required params
		/// Add an option
		template<typename T>
		Option &add_opt(std::string const &n,T &v)
		{
			auto o = std::make_shared<Option>(n);
			o->store(std::make_shared<TValue<T>>(v));
			opt_list.push_back(o);
			return *o;
		}
		/// Perform parsing of the command line 
		void parse(int argc,char *argv[])
		{
			using namespace std;
			last_cmd = argv[0];
			vector<string> pos;
			vector<string> nmd;
			bool name_on = true;
			for (auto o: opt_list) o->_present = false;
			for (int i = 1;i<argc;i++) {
				string s(argv[i]);
				if (name_on && s.substr(0,2)=="--") { // named option
					if (s.length()==2) {
						name_on = false; // turn off name processing
					}
					auto p = s.find('=',2);
					auto id = p==string::npos?s.substr(2):s.substr(2,p-2); // get id of Option
					auto o = find_if( // find name in opt_list
						opt_list.begin(),opt_list.end(),
						[id](std::shared_ptr<Option> n){return n->get_id()==id;}
					);
					if (o==opt_list.end()) throw invalid_argument("unknown option: "+s);
					(*o)->_present = true;
					if ((*o)->vr==Option::VR_NO) {
						if (p!=string::npos) throw invalid_argument(s+" with excessive value");
						(*o)->set_value();
						continue;
					}
					if (p!=string::npos) {
						(*o)->set_value(s.substr(p+1));
						continue;
					}
					if (i+1==argc||string(argv[i+1]).substr(0,2)==string("--")) {
						if ((*o)->vr==Option::VR_YES) throw invalid_argument(s+" missing required value");
						(*o)->set_value();
						continue;
					}
					if ((*o)->vr==Option::VR_YES) {
						(*o)->set_value(argv[i+1]);
						i += 1;
						continue;
					}
					assert((*o)->vr==Option::VR_MAYBE);
					try { // see if next argument can be a value
						(*o)->set_value(argv[i+1]);
						i += 1;
					}
					catch (invalid_argument &e) {
						(*o)->set_value();
					}
					continue;
				}
				if (s[0]=='-') { // keyed option
					throw logic_error("keyed option not implemented");
				}
				pos.push_back(s);
			}
			unsigned pc = 0;
			unsigned sc = 0;
			while (pc<prm_list.size()&&sc<pos.size()) {
				try {
					prm_list[pc]->set_value(pos[sc]);
					sc ++;
				}
				catch (domain_error &e) {
					if (pc<prm_required) throw invalid_argument(pos[sc]+" is invalid for "+prm_list[pc]->get_label());
				}
				pc ++;
			}
			if (pc<prm_required) throw invalid_argument("not enough arguments");
			if (sc<pos.size()) throw invalid_argument("too many arguments");
		}
		/// Check if an Option was present in the last parse()
		/// @param id Name of the Option
		bool opt_is_present(std::string const &id) const
		{
			auto o = find_if( // find name in opt_list
				opt_list.begin(),opt_list.end(),
				[id](std::shared_ptr<Option> opt){return opt->get_id()==id;}
			);
			return o!=opt_list.end() && (*o)->_present;
		}
		/// Get help text
		std::string get_help() const
		{
			using namespace std;
			ostringstream o;
			o << "Usage: " << last_cmd;
			auto rp = prm_required;
			for (auto &p:prm_list) {
				string s = '<'+p->get_label()+'>';
				if (rp) rp--;
				else s = '['+s+']';
				o << ' ' << s;
			}
			o << " [OPTIONS...]\n";
			o << "\n Parameters:\n\n";	
			rp = prm_required;
			for (auto &p:prm_list) {
				string s = "  "+p->get_label();
				s += ' ';
				if (s.length()<16) s.resize(24,' ');
				o << s << p->get_desc();
				if (rp) rp--;
				else o << " [default:" << p->get_str() << ']';
				o << '\n';
			}
			o << "\n Valid options:\n\n";
			for (auto &p:opt_list) {
				string s = "  --"+p->get_id();
				if (p->vr!=Option::VR_NO) {
					string t = '='+p->type_label();
					if (p->vr==Option::VR_MAYBE) t = '['+t+']';
					s += t;
				}
				s += ' ';
				if (s.length()<24) s.resize(24,' ');
				o << s << p->get_desc() << '\n';
			}
			return o.str();
		}
	};
}
