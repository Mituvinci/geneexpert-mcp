#!/usr/bin/env node

import { GoogleGenerativeAI } from '@google/generative-ai';
import dotenv from 'dotenv';

dotenv.config();

const genAI = new GoogleGenerativeAI(process.env.GOOGLE_API_KEY);

async function listModels() {
  console.log('Checking available Gemini models...\n');

  try {
    // Try a simple API call to list models
    const response = await fetch(
      `https://generativelanguage.googleapis.com/v1/models?key=${process.env.GOOGLE_API_KEY}`
    );

    const data = await response.json();

    if (data.models) {
      console.log('Available models:');
      data.models.forEach(model => {
        console.log(`  - ${model.name}`);
        if (model.supportedGenerationMethods?.includes('generateContent')) {
          console.log(`    ✓ Supports generateContent`);
        }
      });
    } else {
      console.log('Error:', data);
    }
  } catch (error) {
    console.error('Error listing models:', error.message);
  }
}

listModels();
