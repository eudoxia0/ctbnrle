/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * see docs/AUTHORS for contributor list
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "streamserial.h"

#include "nullptr03.h"
#include <cstdarg>
#include <iostream>


namespace SERIALNAMESPACE {

    using ::ctbn::nullptr03;

	// helpful formatting fn
	void Indent(std::ostream &os, int indent) {
		for(int i=0;i<indent;i++) os << "\t";
	}

	void IgnoreWS(std::istream &is) {
		while(1) {
			char c = is.peek();
			if (is.fail() || !isspace(c)) return;
			is.get();
		}
	}


	char ReadEscChar(std::istream &is) {
		char c=is.get();
		switch(c) {
			case '\\': return '\\';
			case 'a': return '\a';
			case 'b': return '\b';
			case 'f': return '\f';
			case 'n': return '\n';
			case 'r': return '\r';
			case 't': return '\t';
			case 'v': return '\v';
			default: return c;
		}
	}

	bool ConsumeToken(std::istream &is, const char *tok) {
		int i;
		for(i=0;tok[i]!=0 && tok[i]==is.get();i++)
			;
		if (tok[i]==0) return true;
		for(;i>=0;i--)
			is.unget();
		return false;
	}

	char ReadAmpChar(std::istream &is) {
		if (ConsumeToken(is,"quot;")) return '"';
		if (ConsumeToken(is,"amp;")) return '&';
		if (ConsumeToken(is,"apos;")) return '\'';
		if (ConsumeToken(is,"lt;")) return '<';
		if (ConsumeToken(is,"gt;")) return '>';
		return '&';
	}

	void ReadToken(std::istream &is, std::string &ret) {
		ret.clear();
		IgnoreWS(is);
		bool isquoted = is.peek()=='"';
		if (isquoted) is.get();
		ReadStr(is,ret,isquoted ? "\"" : " \t\n\r\v>\\=");
		if (isquoted) is.get();
	}

	inline bool charin(const char *endchar,char c) {
		return strchr(endchar,c) != nullptr03;
	}

	void ReadStr(std::istream &is, std::string &ret,const char *endchar) {
		ret.clear();
		while(char c=is.peek()) {
			if (charin(endchar,c)) return;
			is.get();
			if (is.fail()) return;
			if (c=='\\') c = ReadEscChar(is);
			else if (c=='&') c = ReadAmpChar(is);
			ret.push_back(c);
		}
	}

	void WriteStr(std::ostream &os, const std::string &s,bool escape) {
		for(std::string::const_iterator i=s.begin();i!=s.end();++i) {
			switch(*i) {
				case '"':
				    os << "&quot;";
				    break;
				case '&':
				    os << "&amp;";
				    break;
				case '\'':
				    os << "&apos;";
				    break;
				case '<':
				    os << "&lt;";
				    break;
				case '>':
				    os << "&gt;";
				    break;
				case '\\':
				    os << "\\\\";
				    break;
				case '\a':
				    (escape ? (os << '\\' << 'a') : (os << '\a'));
					break;
				case '\b':
				    (escape ? (os << '\\' << 'b') : (os << '\b'));
					break;
				case '\f':
				    (escape ? (os << '\\' << 'f') : (os << '\f'));
					break;
				case '\n':
				    (escape ? (os << '\\' << 'n') : (os << '\n'));
					break;
				case '\r':
				    (escape ? (os << '\\' << 'r') : (os << '\r'));
					break;
				case '\t':
				    (escape ? (os << '\\' << 't') : (os << '\t'));
					break;
				case '\v':
				    (escape ? (os << '\\' << 'v') : (os << '\v'));
					break;
				default:
				    os << *i;
				    break;
			}
		}
	}

	void ReadTag(std::istream &is,XMLTagInfo &info) {
		info.attr.clear();
		IgnoreWS(is);
		char c = is.get();
		if (is.fail() || c!='<') throw streamexception("Stream Input Format Error:  expected <");
		if (is.peek()=='\\') {
			info.isstart=false;
			info.isend=true;
			is.get();
		} else {
			info.isstart=true;
			info.isend=false;
		}
		ReadToken(is,info.name);
		IgnoreWS(is);
		if (info.isend) {
			if (is.get()=='>') return;
			throw streamexception("Stream Input Format Error:  expected >");
		}
		while(!is.fail()) {
			char c = is.peek();
			if (c=='\\') {
				is.get();
				info.isend=true;
				if (is.get()=='>') return;
				throw streamexception("Stream Input Format Error:  expected >");
			}
			if (c=='>') {
				is.get();
				return;
			}
			std::string aname;
			ReadToken(is,aname);
			if (aname.empty()) 
				throw streamexception("Stream Input Format Error: tag missing name");
			IgnoreWS(is);
			if (is.peek()=='=') {
				is.get();
				std::string aval;
				ReadToken(is,aval);
				info.attr[aname] = aval;
			} else {
				info.attr[aname] = "1";
			}
			IgnoreWS(is);
		}
		throw streamexception("Stream Input Format Error: unexpected stream end");
	}

	void ReadEndTag(std::istream &is,const char *ename) {
	    XMLTagInfo einfo;
		ReadTag(is,einfo);
		if (!einfo.isend || einfo.name!=ename)
			throw streamexception(std::string("Stream Input Format Error: expected end tag for ")+ename+", received "+(einfo.isend ? "end" : "start")+" tag for "+einfo.name);
	}

	// How to construct a name for a template class
	char *TName(const char *cname, int n, ...) {
		int s = strlen(cname)+1;
		va_list ap;
		va_start(ap,n);
		for(int i=0;i<n;i++)
			s += strlen(va_arg(ap,const char *))+1;
		va_end(ap);
		char *ret = new char[s];
		strcpy(ret,cname);
		va_start(ap,n);
		for(int i=0;i<n;i++) {
			strcat(ret, "."); // slightly inefficient to keep scanning
			strcat(ret, va_arg(ap,const char *)); // for the end, but much
				// more readable
		}
		return ret;
	}

	typedef std::map<std::string,createfntype> allocmaptype;

	allocmaptype &alloctable() {
		static allocmaptype table;
		return table;
	}

	const char *AddAlloc(const char *name, createfntype fn) {
		alloctable()[std::string(name)] = fn;
		return name;
	}

	void *AllocByName(const std::string &name) {
		allocmaptype::iterator i = alloctable().find(name);
		if (i==alloctable().end()) {
			return nullptr03;
		}
		else return i->second();
	}
}
