# SuperDendrix weight function
def SuperW(cases_by_event, profile, samples):
        mut_samples = set(s for cases in cases_by_event for s in cases)
        pos = sum(profile[s] for s in mut_samples if profile[s] > 0)
        neg = sum(abs(profile[s]) for cases in cases_by_event for s in cases)
        return 2 * pos - neg