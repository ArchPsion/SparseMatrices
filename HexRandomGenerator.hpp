#ifndef __HEX_RANDOM_GENERATOR_HPP__
#define __HEX_RANDOM_GENERATOR_HPP__

// Standard Libraries
#include <random>
#include <vector>

// Qt Libraries
#include <QtGlobal>

class HexRandomGenerator
{
	private:
	
		std::random_device				rGen;
		
	public:
	
		inline qint32					getNumberWithinRange(qint32);
		inline quint32					getNumberWithinRange(quint32);
		inline std::vector<qint32>			getNumbersWithinRange(qint32, qint32);
		template<typename Type> inline void		shuffle(std::vector<Type>&);
};

qint32 HexRandomGenerator::getNumberWithinRange(qint32 max)
{
	if (max < 2)
		return 0;
	
	const auto foo = HexRandomGenerator::rGen() % static_cast<quint32>(max);
	return static_cast<qint32>(foo);
}

quint32 HexRandomGenerator::getNumberWithinRange(quint32 max)
{
	if (max < 2u)
		return 0u;
	
	return (HexRandomGenerator::rGen() % max);
}

std::vector<qint32> HexRandomGenerator::getNumbersWithinRange(qint32 quantity, qint32 max)
{	
	if (quantity < 1 or max < 2 or quantity >= max + 1)
		return { };
	
	if (quantity == max)
	{
		auto numbers = std::vector<qint32>(quantity, 0);
		std::iota(numbers.begin(), numbers.end(), 0);
		
		return numbers;
	}
	
	auto numbers = std::vector<qint32>();
	numbers.reserve(quantity*6/5);
	
	if (quantity > max/2) // In this case, you might as well generate the numbers that you WON'T get.
	{
		const auto antiNumbers = HexRandomGenerator::getNumbersWithinRange(max - quantity, max);
		auto index = 0;
		
		for (const auto& val : antiNumbers)
		{
			while (index < val)
			{
				numbers.push_back(index);
				++index;
			}
			
			++index;
		}
		
		while (index < max)
		{
			numbers.push_back(index);
			++index;
		}
		
		return numbers;
	}
	
	auto numbersToProduce = quantity*11/10; // It's unlikely that the n generated numbers will be unique, so we produce a bit more.
	auto uniqueNumbers = 0;
	
	while (numbersToProduce != 0)
	{
		for (auto i = 0; i < numbersToProduce; ++i)
		{
			const auto foo = HexRandomGenerator::getNumberWithinRange(max);
			numbers.emplace_back(foo);
		}
		
		std::sort(numbers.begin(), numbers.end());
		
		auto duplicates = 0;
		auto last = -1;
		
		for (auto& val : numbers)
		{
			if (val == last)
			{
				++duplicates;
				val = -1;
			}
			else
				last = val;
		}
		
		uniqueNumbers = static_cast<qint32>(numbers.size()) - duplicates;
		numbersToProduce = quantity - uniqueNumbers;
		
		while (numbersToProduce < 0)
		{
			auto index = HexRandomGenerator::rGen() % numbers.size(); // We remove this one, or the first non-duplicate element afterward.
			
			while (numbers[index] < 0)
				index = (index + 1u) % numbers.size();
			
			numbers[index] = -1;
			++numbersToProduce;
		}
		
		if (numbersToProduce > 0)
			numbersToProduce = numbersToProduce*11/10;
	}
	
	auto numbersPurged = std::vector<qint32>();
	numbersPurged.reserve(quantity);
	
	for (const auto& val : numbers)
	{
		if (val < 0)
			continue;
		
		numbersPurged.push_back(val);
	}
	
	return numbersPurged;
}

template<typename Type>
void HexRandomGenerator::shuffle(std::vector<Type>& vect)
{
	if (vect.size() > 1u)
		std::shuffle(vect.begin(), vect.end(), HexRandomGenerator::rGen);
}

#endif
