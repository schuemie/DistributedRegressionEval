{DEFAULT @sample_size = 100000}
{DEFAULT @cohort_id = 1}
{DEFAULT @new_cohort_id = 101}

INSERT INTO @cohort_database_schema.@cohort_table (cohort_definition_id,
	subject_id,
	cohort_start_date,
	cohort_end_date)
SELECT @new_cohort_id AS cohort_definition_id,
	subject_id,
	cohort_start_date,
	cohort_end_date
FROM (
	SELECT ROW_NUMBER() OVER (ORDER BY NEWID()) AS rn,
		subject_id,
		cohort_start_date,
		cohort_end_date
	FROM @cohort_database_schema.@cohort_table 
	WHERE cohort_definition_id = @cohort_id
) temp
WHERE rn <= @sample_size;