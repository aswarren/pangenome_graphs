SET ECHO OFF
SET NEWPAGE 0
SET SPACE 0
SET PAGESIZE 0
SET FEEDBACK OFF
SET HEADING OFF
SET HEADSEP OFF
SET TRIMSPOOL ON
SET TAB ON
set sqlprompt ''
set feedback off
set wrap off
set linesize 3000

spool test.txt

select distinct a.ncbi_tax_id, a.rank, tax_path
        from 
        (select ncbi_tax_id, rank, trim(leading ',' from SYS_CONNECT_BY_PATH(ncbi_tax_id,',')) tax_path
        from sres.taxon tx start with tx.ncbi_tax_id in (select ncbi_tax_id from app.genomesummary where genome_info_id in (13058, 28410))
        connect                                                                                                             
        by prior parent_id=taxon_id)  a
        where a.rank = 'genus';

QUIT
