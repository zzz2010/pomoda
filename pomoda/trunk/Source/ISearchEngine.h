#pragma once

class ISearchEngine
{
public:
	ISearchEngine(void);
	virtual ~ISearchEngine(void);
	virtual int searchPattern( char* pattern, int mismatches, vector<VAL>& position)=0;
	virtual VAL getTotalLength()=0;
};
