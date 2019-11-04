//
// Created by meital on 15/11/16.
//

#include <libscapi/include/primitives/Mersenne.hpp>
#include "TemplateField.h"
#include "GF2_8LookupTable.h"


using namespace NTL;

template <>
TemplateField<GF2_8LookupTable>::TemplateField(long fieldParam) {

    elementSizeInBytes = 1;//round up to the next byte
    elementSizeInBits = 8;

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new GF2_8LookupTable(0);
    m_ONE = new GF2_8LookupTable(1);

    GF2_8LookupTable::initTable();
}

template <>
void TemplateField<GF2_8LookupTable>::elementToBytes(unsigned char* elemenetInBytes, GF2_8LookupTable& element){

    memcpy(elemenetInBytes, (byte*)(&element.elem), 1);
}

template <>
GF2_8LookupTable TemplateField<GF2_8LookupTable>::bytesToElement(unsigned char* elemenetInBytes){

    return GF2_8LookupTable((unsigned int)(*(byte*)elemenetInBytes));
}


template <>
GF2_8LookupTable TemplateField<GF2_8LookupTable>::GetElement(long b) {


    if(b == 1)
    {
        return *m_ONE;
    }
    else if(b == 0)
    {
        return *m_ZERO;
    }
    else{
        GF2_8LookupTable element((byte)b);
        return element;
    }
}


template <>
void TemplateField<GF2_8LookupTable>::elementVectorToByteVector(vector<GF2_8LookupTable> &elementVector, vector<byte> &byteVector){

    copy_byte_array_to_byte_vector((byte *)elementVector.data(), elementVector.size()*elementSizeInBytes, byteVector,0);

}

