// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME root_admin
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "ROOT_admin.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_ROOT_admin(void *p = nullptr);
   static void *newArray_ROOT_admin(Long_t size, void *p);
   static void delete_ROOT_admin(void *p);
   static void deleteArray_ROOT_admin(void *p);
   static void destruct_ROOT_admin(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ROOT_admin*)
   {
      ::ROOT_admin *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ROOT_admin >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ROOT_admin", ::ROOT_admin::Class_Version(), "ROOT_admin.h", 36,
                  typeid(::ROOT_admin), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ROOT_admin::Dictionary, isa_proxy, 4,
                  sizeof(::ROOT_admin) );
      instance.SetNew(&new_ROOT_admin);
      instance.SetNewArray(&newArray_ROOT_admin);
      instance.SetDelete(&delete_ROOT_admin);
      instance.SetDeleteArray(&deleteArray_ROOT_admin);
      instance.SetDestructor(&destruct_ROOT_admin);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ROOT_admin*)
   {
      return GenerateInitInstanceLocal((::ROOT_admin*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ROOT_admin*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ROOT_admin::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ROOT_admin::Class_Name()
{
   return "ROOT_admin";
}

//______________________________________________________________________________
const char *ROOT_admin::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ROOT_admin*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ROOT_admin::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ROOT_admin*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ROOT_admin::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ROOT_admin*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ROOT_admin::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ROOT_admin*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ROOT_admin::Streamer(TBuffer &R__b)
{
   // Stream an object of class ROOT_admin.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ROOT_admin::Class(),this);
   } else {
      R__b.WriteClassBuffer(ROOT_admin::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ROOT_admin(void *p) {
      return  p ? new(p) ::ROOT_admin : new ::ROOT_admin;
   }
   static void *newArray_ROOT_admin(Long_t nElements, void *p) {
      return p ? new(p) ::ROOT_admin[nElements] : new ::ROOT_admin[nElements];
   }
   // Wrapper around operator delete
   static void delete_ROOT_admin(void *p) {
      delete ((::ROOT_admin*)p);
   }
   static void deleteArray_ROOT_admin(void *p) {
      delete [] ((::ROOT_admin*)p);
   }
   static void destruct_ROOT_admin(void *p) {
      typedef ::ROOT_admin current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ROOT_admin

namespace {
  void TriggerDictionaryInitialization_libroot_admin_Impl() {
    static const char* headers[] = {
"ROOT_admin.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/farid/softwares/root/include",
"/home/farid/MyRepositories/KineticTheory/kineticTheory/src",
"/usr/local/include",
"/usr/include/hdf5/serial",
"/usr/include",
"/home/farid/softwares/root/include/",
"/home/farid/MyRepositories/KineticTheory/kineticTheory/build/src/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libroot_admin dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$ROOT_admin.h")))  ROOT_admin;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libroot_admin dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "ROOT_admin.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"ROOT_admin", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libroot_admin",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libroot_admin_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libroot_admin_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libroot_admin() {
  TriggerDictionaryInitialization_libroot_admin_Impl();
}
