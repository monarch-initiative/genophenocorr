import hpotk
import typing
import logging

from genophenocorr.model import Phenotype


class PhenotypeValidationException(BaseException):

    def __init__(self, issues):
        self._issues = issues

    @property
    def issues(self):
        return self._issues


class PhenotypeCreator:
    """A class that creates a Phenotype object 

    Methods:
        create_phenotype(term_ids:Iterable[Tuple[str, bool]]): Creates a list of Phenotype objects from a list of tuples.
                                                               Each tuple has the HPO ID and a boolean on if the phenotype is observed.
    """
    def __init__(self, hpo: hpotk.MinimalOntology,
                 validator: hpotk.validate.ValidationRunner):
        """Constructs all necessary attributes for a PhenotypeCreator object

        Args:
            hpo (hpotk.ontology.Ontology): An Ontology object 
            validator (hpotk.validate.ValidationRunner): A ValidationRunner object 
        """
        self._logger = logging.getLogger(__name__)
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._validator = hpotk.util.validate_instance(validator, hpotk.validate.ValidationRunner, 'validator')

    def create_phenotype(self, term_ids: typing.Iterable[typing.Tuple[str, bool]]) -> typing.Tuple[typing.Sequence[Phenotype], typing.Sequence[str]]:
        """Creates a list of Phenotype objects from term IDs and checks if the term IDs satisfy the validation requirements.

        Args:
            term_ids (Iterable[Tuple[str, bool]]): A list of Tuples, structured (HPO IDs, boolean- True if observed)
        Returns:
            A sequence of Phenotype objects
        Error:
            PhenotypeValidationException: An instance of an issue with the ValidationRunner
        """
        terms = []
        ignored_terms = []
        for term_id, observed in term_ids:
            term = self._hpo.get_term(term_id)
            if term is None:
                self._logger.warning("Term %s cannot be found in HPO version %s. It will be ignored.", term_id, self._hpo.version)
                ignored_terms.append(term_id)
            else:
                terms.append((term, observed))
        validation_results = self._validator.validate_all([term[0] for term in terms])
        if validation_results.is_ok:
            return tuple(Phenotype.from_term(term, observed) for term, observed in terms), ignored_terms
        else:
            # We return the messages for now. We may provide more details in future, if necessary.
            issues = [r.message for r in validation_results.results]
            raise PhenotypeValidationException(issues)
