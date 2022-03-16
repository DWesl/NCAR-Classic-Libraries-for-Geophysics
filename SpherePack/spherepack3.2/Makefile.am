ACLOCAL_AMFLAGS=-I m4
AM_LDFLAGS=-no-undefined

lib_LTLIBRARIES=libspherepack.la

libspherepack_la_SOURCES=src/alf.f src/divec.f src/dives.f		\
    src/divgc.f src/divgs.f src/gaqd.f src/geo2math.f src/gradec.f	\
    src/grades.f src/gradgc.f src/gradgs.f src/hrfft.f src/idivec.f	\
    src/idives.f src/idivgc.f src/idivgs.f src/idvtec.f src/idvtes.f	\
    src/idvtgc.f src/idvtgs.f src/igradec.f src/igrades.f		\
    src/igradgc.f src/igradgs.f src/ihgeod.f src/isfvpec.f		\
    src/isfvpes.f src/isfvpgc.f src/isfvpgs.f src/islapec.f		\
    src/islapes.f src/islapgc.f src/islapgs.f src/ivlapec.f		\
    src/ivlapes.f src/ivlapgc.f src/ivlapgs.f src/ivrtec.f		\
    src/ivrtes.f src/ivrtgc.f src/ivrtgs.f src/lfim.f src/lfin.f	\
    src/lfpt.f src/sfvpec.f src/sfvpes.f src/sfvpgc.f src/sfvpgs.f	\
    src/shaec.f src/shaes.f src/shagc.f src/shags.f src/shigc.f		\
    src/shigs.f src/shpe.f src/shpg.f src/shsec.f src/shses.f		\
    src/shsgc.f src/shsgs.f src/slapec.f src/slapes.f src/slapgc.f	\
    src/slapgs.f src/sphcom.f src/sshifte.f src/trssph.f src/trvsph.f	\
    src/vhaec.f src/vhaes.f src/vhagc.f src/vhags.f src/vhsec.f		\
    src/vhses.f src/vhsgc.f src/vhsgs.f src/vlapec.f src/vlapes.f	\
    src/vlapgc.f src/vlapgs.f src/vrtec.f src/vrtes.f src/vrtgc.f	\
    src/vrtgs.f src/vshifte.f src/vtsec.f src/vtses.f src/vtsgc.f	\
    src/vtsgs.f src/math-common.f

check_PROGRAMS=advec helmsph shallow

advec_SOURCES=src/advec.f
advec_LDADD=libspherepack.la

helmsph_SOURCES=src/helmsph.f
helmsph_LDADD=libspherepack.la

shallow_SOURCES=src/shallow.f
shallow_LDADD=libspherepack.la

if HAVE_FFTPACK
libspherepack_la_SOURCES+=src/lfp.f
libspherepack_la_LDFLAGS=-lfftpack
endif

if HAVE_NCAR_GRAPHICS
lib_LTLIBRARIES+=libspherepackvis.la
libspherepackvis_la_SOURCES=src/visequ.f src/visgau.f src/visgeo.f	\
     src/vsurf.f
libspherepackvis_la_LIBADD=libspherepack.la
libspherepackvis_la_LDFLAGS=-lncarg
endif

check_PROGRAMS+=testrvsph tdiv testrssph testsshifte testvshifte testvtsgs \
    tgaqd tgrad tidvt tsha tshpe tshpg tslap tvha tvlap \
    tvrt tvts

TESTS=$(check_PROGRAMS)

testrvsph_SOURCES=test/testrvsph.f
testrvsph_LDADD=libspherepack.la

tdiv_SOURCES=test/tdiv.f
tdiv_LDADD=libspherepack.la

testrssph_SOURCES=test/testrssph.f
testrssph_LDADD=libspherepack.la

testsshifte_SOURCES=test/testsshifte.f
testsshifte_LDADD=libspherepack.la

testvshifte_SOURCES=test/testvshifte.f
testvshifte_LDADD=libspherepack.la

testvtsgs_SOURCES=test/testvtsgs.f
testvtsgs_LDADD=libspherepack.la

tgaqd_SOURCES=test/tgaqd.f
tgaqd_LDADD=libspherepack.la

tgrad_SOURCES=test/tgrad.f
tgrad_LDADD=libspherepack.la

tidvt_SOURCES=test/tidvt.f
tidvt_LDADD=libspherepack.la

tsha_SOURCES=test/tsha.f
tsha_LDADD=libspherepack.la

tshpe_SOURCES=test/tshpe.f
tshpe_LDADD=libspherepack.la

tshpg_SOURCES=test/tshpg.f
tshpg_LDADD=libspherepack.la

tslap_SOURCES=test/tslap.f
tslap_LDADD=libspherepack.la

tvha_SOURCES=test/tvha.f
tvha_LDADD=libspherepack.la

tvlap_SOURCES=test/tvlap.f
tvlap_LDADD=libspherepack.la

tvrt_SOURCES=test/tvrt.f
tvrt_LDADD=libspherepack.la

tvts_SOURCES=test/tvts.f
tvts_LDADD=libspherepack.la
