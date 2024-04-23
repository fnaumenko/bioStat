/**********************************************************
DefRegions.cpp  2023 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 11/21/2023
-------------------------
***********************************************************/

#include "DefRegions.h"
#include "ChromSeq.h"

void DefRegions::Init()
{
	if (IsEmpty())
		if (_cSizes.IsFilled()) {
			// initialize instance from chrom sizes
			if (Chrom::IsSetByUser())
				for (ChromSizes::cIter it = _cSizes.cBegin(); it != _cSizes.cEnd(); it++)
					AddElem(CID(it), Regions(0, _cSizes[CID(it)]));
			else
				AddElem(Chrom::UserCID(), Regions(0, _cSizes[Chrom::UserCID()]));
			//_isEmpty = false;
		}
}

const Regions& DefRegions::operator[] (chrid cID)
{
	if (FindChrom(cID))	return At(cID).Data;
	ChromDefRegions rgns(_cSizes.ServName(cID), _minGapLen);
	if (rgns.Empty())		// file with def regions doesn't exist?
	{
		//_cSizes.IsFilled();
		const string ext = _cSizes.RefExt();
		if (!ext.length())	// no .fa[.gz] file, empty service dir: _cSizes should be initialized by BAM
			return AddElem(cID, Regions(0, _cSizes[cID])).Data;
		//Err(Err::F_NONE, (_cSizes.ServName(cID) + ChromDefRegions::Ext).c_str()).Throw();
		ChromSeq rs(_cSizes.RefName(cID) + ext, rgns, _minGapLen);
	}
	return AddElem(cID, rgns).Data;
}

genlen DefRegions::GenSize() const
{
	genlen gsize = 0;
	for (cIter it = cBegin(); it != cEnd(); it++)
		gsize += Size(it);
	return gsize;
}

chrlen DefRegions::MinSize() const
{
	cIter it = cBegin();
	chrlen	minsize = Size(it);
	for (it++; it != cEnd(); it++)
		if (minsize > Size(it))
			minsize = Size(it);
	return	minsize;
}

#ifdef MY_DEBUG
void DefRegions::Print() const
{
	cout << "DefRegions:\n";
	for (DefRegions::cIter it = cBegin(); it != cEnd(); it++)
		cout << Chrom::TitleName(CID(it))
		<< TAB << Data(it).FirstStart()
		<< TAB << Size(it) << LF;
}
#endif	// MY_DEBUG
